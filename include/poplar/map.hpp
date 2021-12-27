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
#ifndef POPLAR_TRIE_MAP_HPP
#define POPLAR_TRIE_MAP_HPP

#include <array>
#include <iostream>

#include "bit_tools.hpp"
#include "exception.hpp"

namespace poplar {

// This class implements an updatable associative array whose keys are strings.
// The data structure is based on a dynamic path-decomposed trie described in the following paper,
// - "Dynamic Path-Decomposed Tries" available at https://arxiv.org/abs/1906.06015.
template <typename Trie, typename NLM>
class map {
    static_assert(Trie::trie_type_id == NLM::trie_type_id);

  public:
    using this_type = map<Trie, NLM>;
    using trie_type = Trie;
    using value_type = typename NLM::value_type;

    static constexpr auto trie_type_id = Trie::trie_type_id;
    static constexpr uint32_t min_capa_bits = Trie::min_capa_bits;

  public:
    // Generic constructor.
    map() = default;

    // Class constructor. Initially allocates the hash table of length
    // 2**capa_bits.
    explicit map(uint32_t capa_bits, uint64_t lambda = 32) {
        POPLAR_THROW_IF(!is_power2(lambda), "lambda must be a power of 2.");

        is_ready_ = true;
        lambda_ = lambda;
        hash_trie_ = Trie{capa_bits, 8 + bit_tools::ceil_log2(lambda_)};
        label_store_ = NLM{hash_trie_.capa_bits()};
        codes_.fill(UINT8_MAX);
        codes_[0] = static_cast<uint8_t>(num_codes_++);  // terminator
    }

    // Generic destructor.
    ~map() = default;

    // Searches the given key and returns the value pointer if registered;
    // otherwise returns nullptr.
    const value_type* find(const std::string& key) const {
        return find(make_char_range(key));
    }
    const value_type* find(char_range key) const {
        POPLAR_THROW_IF(key.empty(), "key must be a non-empty string.");
        POPLAR_THROW_IF(*(key.end - 1) != '\0', "The last character of key must be the null terminator.");

        if (!is_ready_ or hash_trie_.size() == 0) {
            return nullptr;
        }

        auto node_id = hash_trie_.get_root();

        while (!key.empty()) {
            auto [vptr, match] = label_store_.compare(node_id, key);
            if (vptr != nullptr) {
                return vptr;
            }

            key.begin += match;

            while (lambda_ <= match) {
                node_id = hash_trie_.find_child(node_id, step_symb);
                if (node_id == nil_id) {
                    return nullptr;
                }
                match -= lambda_;
            }

            if (codes_[*key.begin] == UINT8_MAX) {
                // Detecting an useless character
                return nullptr;
            }

            node_id = hash_trie_.find_child(node_id, make_symb_(*key.begin, match));
            if (node_id == nil_id) {
                return nullptr;
            }

            ++key.begin;
        }

        return label_store_.compare(node_id, key).first;
    }

    // Inserts the given key and returns the value pointer.
    value_type* update(const std::string& key) {
        return update(make_char_range(key));
    }
    value_type* update(char_range key) {
        POPLAR_THROW_IF(key.empty(), "key must be a non-empty string.");
        POPLAR_THROW_IF(*(key.end - 1) != '\0', "The last character of key must be the null terminator.");

        if (hash_trie_.size() == 0) {
            if (!is_ready_) {
                *this = this_type{0};
            }
            // The first insertion
            ++size_;
            hash_trie_.add_root();

            if constexpr (trie_type_id == trie_type_ids::FKHASH_TRIE) {
                // assert(hash_trie_.get_root() == label_store_.size());
                return label_store_.append(key);
            }
            if constexpr (trie_type_id == trie_type_ids::BONSAI_TRIE) {
                return label_store_.insert(hash_trie_.get_root(), key);
            }
            // should not come
            assert(false);
        }

        auto node_id = hash_trie_.get_root();

        while (!key.empty()) {
            auto [vptr, match] = label_store_.compare(node_id, key);
            if (vptr != nullptr) {
                return const_cast<value_type*>(vptr);
            }

            key.begin += match;

            while (lambda_ <= match) {
                if (hash_trie_.add_child(node_id, step_symb)) {
                    expand_if_needed_(node_id);
#ifdef POPLAR_EXTRA_STATS
                    ++num_steps_;
#endif
                    if constexpr (trie_type_id == trie_type_ids::FKHASH_TRIE) {
                        assert(node_id == label_store_.size());
                        label_store_.append_dummy();
                    }
                }
                match -= lambda_;
            }

            if (codes_[*key.begin] == UINT8_MAX) {
                // Update table
                restore_codes_[num_codes_] = *key.begin; // 追加
                codes_[*key.begin] = static_cast<uint8_t>(num_codes_++);
                POPLAR_THROW_IF(UINT8_MAX == num_codes_, "");
            }

            if (hash_trie_.add_child(node_id, make_symb_(*key.begin, match))) {
                expand_if_needed_(node_id);
                ++key.begin;
                ++size_;

                if constexpr (trie_type_id == trie_type_ids::FKHASH_TRIE) {
                    assert(node_id == label_store_.size());
                    return label_store_.append(key);
                }
                if constexpr (trie_type_id == trie_type_ids::BONSAI_TRIE) {
                    return label_store_.insert(node_id, key);
                }
                // should not come
                assert(false);
            }

            ++key.begin;
        }

        auto vptr = label_store_.compare(node_id, key).first;
        return vptr ? const_cast<value_type*>(vptr) : nullptr;
    }

    value_type* update_new(char_range key, uint64_t node_id, uint64_t first_match) {
        // std::cout << "--- update_new ---" << std::endl;
        if(!hash_trie_.checK_first_insert()) {
            hash_trie_.set_first_insert(true);
            if constexpr (trie_type_id == trie_type_ids::BONSAI_TRIE) {
                //std::cout << "length : " << key.length() << ", " << key.empty() << std::endl;
                return label_store_.insert_new_table(hash_trie_.get_root(), key); // ルートから出ている
            }
            // should not come
            assert(false);
        }
        
        // 2つ目以降のキー追加の際の処理
        // auto [vptr, match] = label_store_.compare_new_ptrs(node_id, key, first_match);
        // if(vptr == nullptr) return const_cast<value_type*>(vptr);
        bool first_check = true;
        while(!key.empty()) {
            uint64_t match;
            if(first_check) {
                auto [vptr, match_tmp] = label_store_.compare_new_ptrs(node_id, key, first_match);
                if(vptr != nullptr) return const_cast<value_type*>(vptr);
                match = match_tmp;
                first_check = false;
            } else {
                auto [vptr, match_tmp] = label_store_.compare_new_ptrs(node_id, key, 0);
                if(vptr != nullptr) return const_cast<value_type*>(vptr);
                match = match_tmp;
            }
            
            key.begin += match;
            while(lambda_ <= match) {
                if (hash_trie_.add_chid_new_table(node_id, step_symb)) { // step_symbはuint8_tの最大値，つまり255(つまり，ダミーノード)
#ifdef POPLAR_EXTRA_STATS
                    ++num_steps_;
#endif
                }
                match -= lambda_;
            }

            if (codes_[*key.begin] == UINT8_MAX) {
                // Update table
                restore_codes_[num_codes_] = *key.begin;
                codes_[*key.begin] = static_cast<uint8_t>(num_codes_++);
                POPLAR_THROW_IF(UINT8_MAX == num_codes_, "");
            }

            if (hash_trie_.add_chid_new_table(node_id, make_symb_(*key.begin, match))) { // make_symb_は match << 8 | *key.begin
                ++key.begin;

                if constexpr (trie_type_id == trie_type_ids::BONSAI_TRIE) {
                    return label_store_.insert_new_table(node_id, key);
                }
                // should not come
                assert(false);
            }

            ++key.begin;
        }


        auto vptr = label_store_.compare_new_ptrs(node_id, key, 0).first;
        return vptr ? const_cast<value_type*>(vptr) : nullptr;
    }


    // Gets the number of registered keys.
    uint64_t size() const {
        return size_;
    }
    // Gets the capacity of the hash table.
    uint64_t capa_size() const {
        return hash_trie_.capa_size();
    }
#ifdef POPLAR_EXTRA_STATS
    double rate_steps() const {
        return double(num_steps_) / size_;
    }
    uint64_t num_resize() const {
        return hash_trie_.num_resize();
    }
#endif
    uint64_t alloc_bytes() const {
        uint64_t bytes = 0;
        bytes += hash_trie_.alloc_bytes();
        bytes += label_store_.alloc_bytes();
        bytes += codes_.size();
        return bytes;
    }

    // 追加
    // 特定のノードから、get_root()までの文字列を復元する
    std::string restore_insert_string(uint64_t node_id) {
        std::string insert_string = ""; // ここに文字列を格納して、新しい辞書に挿入する

        // とりあえず、文字列を復元する
        auto fs = label_store_.return_string(node_id);
        if(fs == nullptr) return insert_string;
        for(uint64_t i=0;; i++) {
            if(fs[i] == 0x00) break;
            insert_string += fs[i];
        }

        // get_root()まで、文字列を復元する
        while(node_id != hash_trie_.get_root()) {
            auto [parent, symb] = hash_trie_.get_parent_and_symb(node_id); // 親ノードとsymbを取得
            auto [c, match] = restore_symb_(symb); // symbから、遷移に失敗した箇所とlabelを取得する

            insert_string = c + insert_string;

            uint64_t dummy_step = 0; // ダミーノードの数を数える
            while(1) {
                fs = label_store_.return_string(parent);
                if(fs == nullptr) {
                    dummy_step++;
                    auto [tmp1, tmp2] = hash_trie_.get_parent_and_symb(parent);
                    parent = tmp1;
                } else {
                    break;
                }
            }

            match += dummy_step * lambda_; // スキップした回数分足してあげる

            if(match != 0) {
                fs = label_store_.return_string(parent);
                std::string tmp_str = "";
                for(uint64_t j=0; j < match; j++) {
                    tmp_str += fs[j];
                }
                insert_string = tmp_str + insert_string;
            }

            node_id = parent;
        }
        return insert_string;
    }

    // 追加
    // cp_orderから新しいdynpdtを構築する
    void restructure_centroid_path_order(const std::vector<uint64_t>& cp_order, const std::vector<bool>& check_bottom) {
        // 情報を格納するための新しい辞書を用意する
        hash_trie_.expand_new_table();
        label_store_.expand_ptrs();

        for(auto node_id : cp_order) {
            if(check_bottom[node_id]) { // ノード上でみた際に、一番そこの部分(起点)
                std::string insert_key = restore_insert_string(node_id);
                int* ptr = update_new(make_char_range(insert_key), hash_trie_.get_root(), 0);
                *ptr = 1;
            } else { // 起点にしたノードと比較することで、格納場所を求める
                auto fs = label_store_.return_string(node_id);
                if(fs == nullptr) return;
                std::string insert_key = restore_insert_string(node_id);
                int* ptr = update_new(make_char_range(insert_key), hash_trie_.get_root(), 0);
                *ptr = 1;
            }
        }

        // 最後にmoveさせて、終了
        hash_trie_.move_table();
        label_store_.move_ptrs();
        hash_trie_.set_first_insert(false);
    }

    // 追加
    // plain_bonsai_trieないから、calc_topoを呼び出す
    void call_topo() {
        auto [cp_order, check_bottom] = hash_trie_.calc_topo(restore_codes_);

        // 新しい辞書を作成し、登録する
        restructure_centroid_path_order(cp_order, check_bottom);
    }

    void show_stats(std::ostream& os, int n = 0) const {
        auto indent = get_indent(n);
        show_stat(os, indent, "name", "map");
        show_stat(os, indent, "lambda", lambda_);
        show_stat(os, indent, "size", size());
        show_stat(os, indent, "alloc_bytes", alloc_bytes());
#ifdef POPLAR_EXTRA_STATS
        show_stat(os, indent, "rate_steps", rate_steps());
#endif
        show_member(os, indent, "hash_trie_");
        hash_trie_.show_stats(os, n + 1);
        show_member(os, indent, "label_store_");
        label_store_.show_stats(os, n + 1);
    }

    map(const map&) = delete;
    map& operator=(const map&) = delete;

    map(map&&) noexcept = default;
    map& operator=(map&&) noexcept = default;

  private:
    static constexpr uint64_t nil_id = Trie::nil_id;
    static constexpr uint64_t step_symb = UINT8_MAX;  // (UINT8_MAX, 0)

    bool is_ready_ = false;
    uint64_t lambda_ = 32;

    Trie hash_trie_;
    NLM label_store_;
    std::array<uint8_t, 256> codes_ = {};
    std::array<uint8_t, 256> restore_codes_ = {}; // 追加
    uint32_t num_codes_ = 0;
    uint64_t size_ = 0;
#ifdef POPLAR_EXTRA_STATS
    uint64_t num_steps_ = 0;
#endif

    uint64_t make_symb_(uint8_t c, uint64_t match) const {
        assert(codes_[c] != UINT8_MAX);
        return static_cast<uint64_t>(codes_[c]) | (match << 8);
    }

    // 追加
    // labelからcとmatchを復元する
    std::pair<char, uint64_t> restore_symb_(uint64_t label) const  {
        return std::pair{char(restore_codes_[label % 256]), label/256};
    }

    void expand_if_needed_(uint64_t& node_id) {
        if constexpr (trie_type_id == trie_type_ids::BONSAI_TRIE) {
            if (!hash_trie_.needs_to_expand()) {
                return;
            }
            auto node_map = hash_trie_.expand();
            node_id = node_map[node_id];
            label_store_.expand(node_map);
        }
    }
};

}  // namespace poplar

#endif  // POPLAR_TRIE_MAP_HPP
