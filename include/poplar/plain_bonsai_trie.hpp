/**
 * MIT License
 *
 * Copyright (c) 2018–2019 Shunsuke Kanda
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef POPLAR_TRIE_PLAIN_BONSAI_TRIE_HPP
#define POPLAR_TRIE_PLAIN_BONSAI_TRIE_HPP

#include "bit_tools.hpp"
#include "bit_vector.hpp"
#include "compact_vector.hpp"
#include "hash.hpp"

#include <queue>
#include <algorithm>

namespace poplar {

template <uint32_t MaxFactor = 90, typename Hasher = hash::vigna_hasher>
class plain_bonsai_trie {
  private:
    static_assert(0 < MaxFactor and MaxFactor < 100);

  public:
    static constexpr uint64_t nil_id = UINT64_MAX;
    static constexpr uint32_t min_capa_bits = 16;

    static constexpr auto trie_type_id = trie_type_ids::BONSAI_TRIE;

  public:
    plain_bonsai_trie() = default;

    plain_bonsai_trie(uint32_t capa_bits, uint32_t symb_bits) {
        capa_size_ = size_p2{std::max(min_capa_bits, capa_bits)};
        symb_size_ = size_p2{symb_bits};
        max_size_ = static_cast<uint64_t>(capa_size_.size() * MaxFactor / 100.0);
        table_ = compact_vector{capa_size_.size(), capa_size_.bits() + symb_size_.bits()};
    }

    ~plain_bonsai_trie() = default;

    uint64_t get_root() const {
        assert(size_ != 0);
        return 1;
    }

    void add_root() {
        assert(size_ == 0);
        size_ = 1;
    }

    uint64_t find_child(uint64_t node_id, uint64_t symb) const {
        assert(node_id < capa_size_.size());
        assert(symb < symb_size_.size());

        if (size_ == 0) {
            return nil_id;
        }

        uint64_t key = make_key_(node_id, symb);
        assert(key != 0);

        for (uint64_t i = Hasher::hash(key) & capa_size_.mask();; i = right_(i)) {
            if (i == 0) {
                // table_[0] is always empty so that table_[i] = 0 indicates to be empty.
                continue;
            }
            if (i == get_root()) {
                continue;
            }
            if (table_[i] == 0) {
                // encounter an empty slot
                return nil_id;
            }
            if (table_[i] == key) {
                return i;
            }
        }
    }

    bool add_child(uint64_t& node_id, uint64_t symb) {
        assert(node_id < capa_size_.size());
        assert(symb < symb_size_.size());

        uint64_t key = make_key_(node_id, symb);
        assert(key != 0);

        for (uint64_t i = Hasher::hash(key) & capa_size_.mask();; i = right_(i)) {
            if (i == 0) {
                // table_[0] is always empty so that any table_[i] = 0 indicates to be empty.
                continue;
            }

            if (i == get_root()) {
                continue;
            }

            if (table_[i] == 0) {
                // this slot is empty
                if (size_ == max_size_) {
                    return false;  // needs to expand
                }

                table_.set(i, key);

                ++size_;
                node_id = i;

                return true;
            }

            if (table_[i] == key) {
                node_id = i;
                return false;  // already stored
            }
        }
    }

    std::pair<uint64_t, uint64_t> get_parent_and_symb(uint64_t node_id) const {
        assert(node_id < capa_size_.size());

        uint64_t key = table_[node_id];
        if (key == 0) {
            // root or not exist
            return {nil_id, 0};
        }
        // Returns pair (parent, label)
        return std::make_pair(key >> symb_size_.bits(), key & symb_size_.mask());
    };

    class node_map {
      public:
        node_map() = default;

        node_map(compact_vector&& map, bit_vector&& done_flags)
            : map_{std::move(map)}, done_flags_{std::move(done_flags)} {}

        ~node_map() = default;

        uint64_t operator[](uint64_t i) const {
            return done_flags_[i] ? map_[i] : UINT64_MAX;
        }

        uint64_t size() const {
            return map_.size();
        }

        node_map(const node_map&) = delete;
        node_map& operator=(const node_map&) = delete;

        node_map(node_map&& rhs) noexcept = default;
        node_map& operator=(node_map&& rhs) noexcept = default;

      private:
        compact_vector map_;
        bit_vector done_flags_;
    };

    bool needs_to_expand() const {
        return max_size() <= size();
    }

    node_map expand() {
        plain_bonsai_trie new_ht{capa_bits() + 1, symb_size_.bits()};
        new_ht.add_root();

#ifdef POPLAR_EXTRA_STATS
        new_ht.num_resize_ = num_resize_ + 1;
#endif

        bit_vector done_flags(capa_size());
        done_flags.set(get_root());

        table_.set(get_root(), new_ht.get_root());

        std::vector<std::pair<uint64_t, uint64_t>> path;
        path.reserve(256);

        // 0 is empty, 1 is root
        for (uint64_t i = 2; i < table_.size(); ++i) {
            if (done_flags[i] || table_[i] == 0) {
                // skip already processed or empty elements
                continue;
            }

            path.clear();
            uint64_t node_id = i;

            do {
                auto [parent, label] = get_parent_and_symb(node_id);
                assert(parent != nil_id);
                path.emplace_back(std::make_pair(node_id, label));
                node_id = parent;
            } while (!done_flags[node_id]);

            uint64_t new_node_id = table_[node_id];

            for (auto rit = std::rbegin(path); rit != std::rend(path); ++rit) {
                new_ht.add_child(new_node_id, rit->second);
                table_.set(rit->first, new_node_id);
                done_flags.set(rit->first);
            }
        }

        node_map node_map{std::move(table_), std::move(done_flags)};
        std::swap(*this, new_ht);

        return node_map;
    }

    // 追加
    // トポロジカルソートを求めるときに情報を記憶しておく型
    struct info_fp {
        uint64_t match; // 分岐位置
        uint64_t cnt; // 分岐位置以降の葉の数
        std::vector<uint64_t> children; // 子の集合

        info_fp() : match(0), cnt(0) {}

        info_fp(uint64_t m, uint64_t c) : match(m), cnt(c) {}
    };

    // cpを求めるための関数
    template <class CP, class CL, class FP>
    void find_centrpod_path(CP& cp, const CL& cl, FP& fp, uint64_t node_id) {
        if(fp[node_id].size() == 0) { // 一番まで、たどり着いた時の処理
            // cp.push_back(node_id);
            return;
        }

        // lambda関数の定義(shelterの中身を求めるための関数)
        // auto get_shelter = [&](uint64_t pos) { // posはfp[][pos]に対応
        //     std::vector<std::pair<uint64_t, uint64_t>> shelter;
        //     for(uint64_t i=0; i < fp[node_id][pos].children.size(); i++) {
        //         uint64_t next_node_id = fp[node_id][pos].children[i];
        //         shelter.push_back({next_node_id, cl[next_node_id]});
        //     }
        //     std::sort(shelter.begin(), shelter.end(), [] (auto l, auto r) { // 葉が多い順にソート
        //         return l.second > r.second;
        //     });
        //     return shelter;
        // };

        bool check_zero = false;
        if(fp[node_id][0].match == 0) check_zero = true;

        // 0分岐とそれ以外の個数を比べて，多い方から順に処理する
        if(check_zero) { // 0分岐がある時
            uint64_t sum = 0;
            std::sort(fp[node_id].begin()+1, fp[node_id].end(), [] (auto l, auto r) {
                return l.cnt > r.cnt;
            });
            for(uint64_t i=1; i < fp[node_id].size(); i++) {
                sum += fp[node_id][i].cnt;
            }
            if(fp[node_id][0].cnt > sum) { // 最初に0分岐を処理する
                // それぞれの子を持ってくる、その中から、cntが多い順に処理する
                std::vector<std::pair<uint64_t, uint64_t>> shelter;
                for(uint64_t i=0; i < fp[node_id][0].children.size(); i++) {
                    uint64_t next_node_id = fp[node_id][0].children[i];
                    shelter.push_back({next_node_id, cl[next_node_id]});
                }
                std::sort(shelter.begin(), shelter.end(), [] (auto l, auto r) {
                    return l.second > r.second;
                });
                // auto shelter = get_shelter(0);
                for(auto s : shelter) {
                    find_centrpod_path(cp, cl, fp, s.first);
                    cp.push_back(s.first);
                }

                // 1分岐以上の処理
                for(uint64_t i=1; i < fp[node_id].size(); i++) {
                    shelter.clear();
                    for(uint64_t j=0; j < fp[node_id][i].children.size(); j++) {
                        uint64_t next_node_id = fp[node_id][i].children[j];
                        shelter.push_back({next_node_id, cl[next_node_id]});
                    }
                    std::sort(shelter.begin(), shelter.end(), [] (auto l, auto r) {
                        return l.second > r.second;
                    });
                    // shelter = get_shelter(i);
                    for(auto s : shelter) {
                        find_centrpod_path(cp, cl, fp, s.first);
                        cp.push_back(s.first);
                    }
                }
            } else { // 最後に0分岐を処理する
                // 1分岐以上の処理
                std::vector<std::pair<uint64_t, uint64_t>> shelter;
                for(uint64_t i=1; i < fp[node_id].size(); i++) {
                    shelter.clear();
                    for(uint64_t j=0; j < fp[node_id][i].children.size(); j++) {
                        uint64_t next_node_id = fp[node_id][i].children[j];
                        shelter.push_back({next_node_id, cl[next_node_id]});
                    }
                    std::sort(shelter.begin(), shelter.end(), [] (auto l, auto r) {
                        return l.second > r.second;
                    });
                    // auto shelter = get_shelter(i);
                    for(auto s : shelter) {
                        find_centrpod_path(cp, cl, fp, s.first);
                        cp.push_back(s.first);
                    }
                }

                // 0分岐の処理
                shelter.clear();
                for(uint64_t i=0; i < fp[node_id][0].children.size(); i++) {
                    uint64_t next_node_id = fp[node_id][0].children[i];
                    shelter.push_back({next_node_id, cl[next_node_id]});
                }
                std::sort(shelter.begin(), shelter.end(), [] (auto l, auto r) {
                    return l.second > r.second;
                });
                // shelter = get_shelter(0);
                for(auto s : shelter) {
                    find_centrpod_path(cp, cl, fp, s.first);
                    cp.push_back(s.first);
                }
            }
        } else { // 0分岐がないとき
            std::sort(fp[node_id].begin(), fp[node_id].end(), [] (auto l, auto r) {
                return l.cnt > r.cnt;
            });
            std::vector<std::pair<uint64_t, uint64_t>> shelter;
            for(uint64_t i=0; i < fp[node_id].size(); i++) {
                shelter.clear();
                for(uint64_t j=0; j < fp[node_id][i].children.size(); j++) {
                    uint64_t next_node_id = fp[node_id][i].children[j];
                    shelter.push_back({next_node_id, cl[next_node_id]});
                }
                std::sort(shelter.begin(), shelter.end(), [] (auto l, auto r) {
                    return l.second > r.second;
                });
                // auto shelter = get_shelter(i);
                for(auto s : shelter) {
                    find_centrpod_path(cp, cl, fp, s.first);
                    cp.push_back(s.first);
                }
            }
        }
        if(node_id == get_root()) cp.push_back(node_id);
    }

    // 追加
    // トポロジカルソートで求める
    template <typename C>
    void calc_topo(const C& restore_codes_) {
        std::cout << "--- calc_topo ---" << std::endl;
        uint64_t table_size = table_.size();
        std::cout << "table_size : " << table_size << std::endl;
        std::vector<std::pair<uint64_t, uint64_t>> parent(table_size); // 親の位置を保存するための配列(位置、分岐位置)
        std::vector<uint64_t> partial_num(table_size, 0); // 子の数を格納するための配列(自身の数も含む)
        
        std::vector<std::vector<info_fp>> fork_pos(table_size);// std::vector<std::vector<std::pair<uint64_t, uint64_t>>> fork_pos(table_size); // それぞれの分岐位置で個数を求めるためのもの(分岐位置)
        std::vector<uint64_t> cnt_leaf(table_size, 0); // それぞれのノードから繋がっている葉ノードの数をカウント

        // O(n)で、親の位置、子の数(ノード番号も)を数える
        for(uint64_t i=0; i < table_size; i++) {
            if(table_[i] != 0) {
                auto [p, label] = get_parent_and_symb(i); // 親と遷移情報の取得
                auto [c, match] = std::pair{uint8_t(restore_codes_[label % 256]), label/256}; // 遷移文字と分岐位置を取得
                partial_num[i] += 1; // 自身の数をカウントする
                partial_num[p] += 1; // 子から親の数をカウントする
                parent[i].first = p;
                parent[i].second = match;
            }
        }

        // queueを使用して、一番下のものから処理していく(CPを求める)
        std::queue<uint64_t> que;
        // 対象のデータを集めてくる
        for(uint64_t i=0; i < table_size; i++) {
            if(partial_num[i] == 1) que.push(i);
        }
        // 数を数える
        uint64_t loop_cnt = 0;
        while(!que.empty()) { // 追加された順に処理
            uint64_t q = que.front();
            que.pop();
            cnt_leaf[q] += 1;
            auto [node_id, match] = parent[q];
            uint64_t p = node_id;
            cnt_leaf[p] += cnt_leaf[q];
            // fork_posに対して、分岐の位置に対して、葉がいくつあるのかをカウント
            bool flag = false;
            for(uint64_t j=0; j < fork_pos[p].size(); j++) {
                if(fork_pos[p][j].match == match) {
                    flag = true;
                    fork_pos[p][j].cnt += cnt_leaf[q];
                    fork_pos[p][j].children.push_back(q);
                    break;
                } else if(fork_pos[p][j].match > match) {
                    break;
                }
            }
            if(!flag) {
                uint64_t pos = fork_pos[p].size();
                fork_pos[p].push_back(info_fp{match, cnt_leaf[q]});
                fork_pos[p][pos].children.push_back(q);
                while(1) {
                    if(pos == 0) break;
                    if(fork_pos[p][pos].match < fork_pos[p][pos-1].match) {
                        // 入れ替え
                        auto tmp = fork_pos[p][pos];
                        fork_pos[p][pos] = fork_pos[p][pos-1];
                        fork_pos[p][pos-1] = tmp;
                        pos--;
                    } else {
                        break;
                    }
                }
            }

            partial_num[p]--;
            if(partial_num[p] == 1) { // 子供の処理がすべて終了すると追加
                que.push(p);
            }
            loop_cnt++;
        }
        std::cout << "loop_cnt : " << loop_cnt << std::endl;
        std::cout << "cnt_leaf(get_root) : " << cnt_leaf[get_root()] << std::endl;

        // Centroid Pathを求める
        // とりあえず大きい順に配列に格納してみる
        std::vector<uint64_t> cp_order;
        find_centrpod_path(cp_order, cnt_leaf, fork_pos, get_root());
        std::cout << "cp_order_size : " << cp_order.size() << std::endl;

        // uint64_t cnt = 0;
        // for(uint64_t i=0; i < table_size; i++) {
        //     if(parent[i].first != 0) cnt++;
        // }
        // std::cout << "cnt : " << cnt << std::endl;

        // std::cout << "root_num : " << fork_pos[get_root()].size() << std::endl;
        // for(auto p : fork_pos[get_root()]) {
        //     std::cout << p.match << ", " << p.cnt;
        //     // for(auto c : p.children) {
        //     //     std::cout << ", " << c;
        //     // }
        //     std::cout << std::endl;
        // }
    }

    // # of registerd nodes
    uint64_t size() const {
        return size_;
    }
    uint64_t max_size() const {
        return max_size_;
    }
    uint64_t capa_size() const {
        return capa_size_.size();
    }
    uint32_t capa_bits() const {
        return capa_size_.bits();
    }
    uint64_t symb_size() const {
        return symb_size_.size();
    }
    uint32_t symb_bits() const {
        return symb_size_.bits();
    }
#ifdef POPLAR_EXTRA_STATS
    uint64_t num_resize() const {
        return num_resize_;
    }
#endif
    uint64_t alloc_bytes() const {
        return table_.alloc_bytes();
    }

    void show_stats(std::ostream& os, int n = 0) const {
        auto indent = get_indent(n);
        show_stat(os, indent, "name", "plain_bonsai_trie");
        show_stat(os, indent, "factor", double(size()) / capa_size() * 100);
        show_stat(os, indent, "max_factor", MaxFactor);
        show_stat(os, indent, "size", size());
        show_stat(os, indent, "alloc_bytes", alloc_bytes());
        show_stat(os, indent, "capa_bits", capa_bits());
        show_stat(os, indent, "symb_bits", symb_bits());
#ifdef POPLAR_EXTRA_STATS
        show_stat(os, indent, "num_resize", num_resize_);
#endif
    }

    plain_bonsai_trie(const plain_bonsai_trie&) = delete;
    plain_bonsai_trie& operator=(const plain_bonsai_trie&) = delete;

    plain_bonsai_trie(plain_bonsai_trie&&) noexcept = default;
    plain_bonsai_trie& operator=(plain_bonsai_trie&&) noexcept = default;

  private:
    compact_vector table_;
    uint64_t size_ = 0;  // # of registered nodes
    uint64_t max_size_ = 0;  // MaxFactor% of the capacity
    size_p2 capa_size_;
    size_p2 symb_size_;
#ifdef POPLAR_EXTRA_STATS
    uint64_t num_resize_ = 0;
#endif

    uint64_t make_key_(uint64_t node_id, uint64_t symb) const {
        return (node_id << symb_size_.bits()) | symb;
    }
    uint64_t right_(uint64_t slot_id) const {
        return (slot_id + 1) & capa_size_.mask();
    }
};

}  // namespace poplar

#endif  // POPLAR_TRIE_PLAIN_BONSAI_TRIE_HPP
