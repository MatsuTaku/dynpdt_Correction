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
#include <iostream>
#include <map>
#include <tuple>
#include <iostream>
#include <chrono>
#include <fstream>
#include <random>
#include <poplar.hpp>
#include <algorithm>
#include <unistd.h>

namespace {

using value_type = int;

// 比較した際の深度を表示
void show_hist(const std::map<int, int>& mp) {
    for(auto p : mp) {
        std::cout << p.first << ", " << p.second << std::endl;
    }
}

class Stopwatch {
  using hrc = std::chrono::high_resolution_clock;
  hrc::time_point start_;
 public:
  Stopwatch() : start_(hrc::now()) {}
  auto time_process() const {
    return hrc::now() - start_;
  }
  double get_sec() const {
    return std::chrono::duration<double>(time_process()).count();
  }
  double get_milli_sec() const {
    return std::chrono::duration<double, std::milli>(time_process()).count();
  }
  double get_micro_sec() const {
    return std::chrono::duration<double, std::micro>(time_process()).count();
  }
};

template <class Process>
double milli_sec_in(Process process) {
  Stopwatch sw;
  process();
  return sw.get_milli_sec();
}

inline uint64_t get_process_size() {
    FILE* fp = std::fopen("/proc/self/statm", "r");
    uint64_t dummy(0), vm(0);
    std::fscanf(fp, "%ld %ld ", &dummy, &vm);  // get resident (see procfs)
    std::fclose(fp);
    return vm * ::getpagesize();
}

std::vector<std::string> keys;
std::vector<std::string> time_keysets;

// キー検索速度を計測するためのkeysetsの作成
void MakeTimeKeysets(uint64_t size) {
    std::random_device rnd;
    time_keysets.clear();
    for(int i=0; i < 100000; i++) {
        int v = rnd() % size;
        time_keysets.push_back((keys[v]));
    }
}

template <typename MAP>
double AveFind(const MAP &dyn_, uint64_t size) {
    double time_sum = 0.0;
    for(int i=0; i < 10; i++) {
        MakeTimeKeysets(size);
        //std::sort(time_keysets.begin(), time_keysets.end());
        Stopwatch sw;
        for (int j = 0; j < 100000; j++) {
            int cnt = 0;
            const int *ptr = dyn_.find(time_keysets[j]);
            if (not (ptr != nullptr and *ptr == 1)) {
                std::cerr << "Failed to find " << time_keysets[j] << std::endl;
                return -1;
            }
        }
        time_sum += sw.get_milli_sec();
    }
    time_sum /= 10.0;
    return time_sum;
}

void FileRead() {
    std::string input_name = "../../../dataset/Titles-enwiki.txt";
    // std::string input_name = "../../../dataset/DS5";
    // std::string input_name = "../../../dataset/GeoNames.txt";
    // std::string input_name = "../../../dataset/AOL.txt";

    // std::string input_name = "../../../dataset/enwiki-20150205.line";
    // std::string input_name = "../../../dataset/wordnet-3.0-word";
    
    std::ifstream ifs(input_name);
    if (!ifs) {
        std::cerr << "File not found input file: "<< input_name << std::endl;
        exit(0);
    }
    std::cout << "dataset : " << input_name.substr(17) << std::endl;
    for (std::string s; std::getline(ifs, s); ) {
        keys.push_back(s);
    }
}

template<class Map>
void test()
{
    FileRead();

    uint32_t capa_bits = 16;
    uint64_t lambda = 32;
    const auto num_keys = static_cast<int>(keys.size());
    std::cout << "num_keys : " << num_keys << std::endl;

    auto map = std::make_unique<Map>(capa_bits, lambda);

    for(int i=0; i < num_keys; i++) {
        *map->update(keys[i]) = 1;
        for(int j=0; j <=i; j++) {
            int cnt = 0;
            auto ptr = map->find(keys[j], cnt);
            if(not (ptr != nullptr and *ptr == 1)) {
                std::cout << "failed : " << keys[j] << std::endl;
                std::cout << "loop_cnt : " << i << std::endl;
                std::cout << "ptr : " << *ptr << std::endl;
                return;
            }
        }
    }
    std::cout << "ok." << std::endl;
}

} // namespace

int main() {
    FileRead();
    const auto num_keys = static_cast<int>(keys.size());
    std::cout << "num_keys : " << num_keys << std::endl;

    auto begin_size = get_process_size();

    poplar::plain_bonsai_map<int> map;

    Stopwatch sw;
    for(int i=0; i < num_keys; i++) {
        int* ptr = map.update(keys[i]);
        *ptr = 1;
    }

    auto time = sw.get_milli_sec();
    std::cout << "----------------------" << std::endl;
    bool test_check = true;
    for(int i=0; i < num_keys; i++) {
        const int *ptr = map.find(keys[i]);
        if(not (ptr != nullptr and *ptr == 1)) {
            std::cout << "search_failed : " << keys[i] << std::endl;
            std::cout << "ptr : " << *ptr << std::endl;
            test_check = false;
            return 0;
        }
    }
    std::cout << (test_check ? "ok." : "failed.") << std::endl;
    
    auto ram_size = get_process_size() - begin_size;
    // auto search_time = AveFind(map, keys.size());
    std::cout << "Build time(m): " << time / 1000.0 << std::endl;
    std::cout << "ram_size : " << ram_size << std::endl;
    // std::cout << "time_search : " << search_time << std::endl;
    std::cout << "capa_size : " << map.capa_size() << std::endl;

    return 0;
}
