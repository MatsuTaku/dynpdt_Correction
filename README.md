# 添削用のレポジトリ

## コードの動かし方
> git clone https://github.com/kawanishik/dynpdt_Correction.git<br>
> cd dynpdt_Correction<br>
> cd build<br>
> cmake ..<br>
> make<br>
> cd sample<br>
> ./sample

## 変更箇所について
- 基本的に変更箇所には[// 追加]をコメントしています
1. sample/sample.cpp
   - キーワード辞書に単語を登録したり、検索したりしています
   - 172行目のmap.call_topo();で，該当の関数を呼び出しています

2. include/poplar/map.cpp
  - 222行目で、plain_bonsai_trie.hppから該当関数を呼び出している

3. include/poplar/plain_bonsai_trie.hpp
   - このファイルで、Centroid Pathを求めています
   - 関数calc_topo[ 242行目 ]
     - ノード間のつながりを求める関数
     - ノードのどの位置から、葉が何個あるのかを研鑽している
   - 関数find_centrpod_path[ 367行目 ]
     - 求めた繋がりから、一番深い順から、CPとして、求める