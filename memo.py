


### 20211010_sawada
### 一つ目の反射が本来存在していなさそうなのに存在している done
### num_ch = 30 以上としてしまうと謎の発散現象が発生する done てかしょうがない
### array[j1][i1] < func_fm(array[j1][i1]) が全てになってしまう done
### 「10_purdue_RDE_fresh_jouge.py」がそもそも上下対象にならない done

### 20211011_sawada
### 「rde_l よりも上」の点を消す
### 「os よりも下」の点を消す -> でも消すべきかどうかは要検討
### 格子を描く
### os の後の計算

### 20211019_sawada
### 昔の計算結果を引張ってきた
### mach_number & pressure 分かればその他の値を特定できるつもり
### y+ & y- は更新しているがあくまで特性線の横断する点の座標でしかない，initial_point の座標が変化するわけではない
### 勘違いしてはならないことは「点の値」と「特性線上の値」は分けて考えないといけないということ
### 一度 example 17.1 の例題を「example_17dot1.py」で再現してみて done

### 20211021_sawada
### 対象の要素のarray_alpha[j][i]（他にもあるかもしれないけれども）計算してない
### だからこそ2つ目の要素で計算が進まない
### でも，このまま計算がうまくいくのであれば，x = 1.6, y = 0.58 が確定してしまうの？　どこかでリセットされてる？
### どうして長さの次元を変えても同じなのか？
### x[-], y[-] で本当にいいのかどうか？　x[m], y[m] では？ 
### 「example_17dot1.py」は x[m], y[m] -> x*0.1 [m], y*0.1 [m] としたら *0.1 の座標の結果が返ってくる


### 20211022_sawada
### 2列目で計算がうまくいかないのは列方向に更新していない何かがあるのでは？

### 20211023_sawada 
### predictor の state_4 のパラメーターが全部大きい -> 違う -> 正しい，S_add の補正のせいで勘違いした done
### corrector の p_plus, theta_plus, V_plus 等を使用し始めると state 4 の値が大きくなる -> 違う done
### R : 気体定数，R_o : rho*V のわけのわからない係数


### 20211026_sawada
### 最初の流線の角度のΔthetaめちゃくちゃ小さくしてもそんなに小さくならない
### おそらくこれ以上はデトネーション波の条件次第な気がする
### 流速自体，この後の計算で小さくなるわけやし，最初の点がもっとデトネーション波に近づく？ 近づかない
### 普通の等間隔でとりあえず進めてほしい

### 20211029_sawada
### CJ_speed なんて今のところ，未燃混合気相の部分でしか関わってこない
### -> 最初の一点目がデトネーション波に近づくなんてことはない
### upper free boundary condition の array[0][i] にて誤って S_add を加えてしまっていた
### T_plus の計算の (array_x[0][i]-array_x[1][i-1]) が大きい -> T_plus が小さくなってしまう？
### でも S_plus も変に小さい
### あんまりうまく最初の角度の変化量を大きくできないけども
### おそらく最初のいくつかの流れ線の角度の変化を小さくして，それ以降等差の角度変化を与えている気がするFievson

### 20211031_sawada
### 境界に近い要素のlambdaって何に依存して変化するんでしょうか？
### デトネーション波付近の初期値を計算する際に用いるfunc_neu2Mの許容誤差がほしい初期値と同じオーダーにならないように
### 垂直の特性線になる場合の対応策が必要
### func_neu2M でのepsが小さいほど，M=1のときのMがそもそも精度が悪くなり（やっぱり関係ない），精度の悪い傾きが発生
### 原点がどこであれ傾きが無限になる

### 20211101_sawada
### そもそも今更ながら上下対象になっていない

### 20211102_sawada
### T+, T- のところ保存されているんだけどなー
### やられた...
### 関数 S の中の delta はやはり二次元平面分布(delta=0)，二次元錐状分布(delta=1)
### 澤田これまで delta = 1 で計算してしまった
### delta = 0 にしたらしっかりと上下対称の特性線を計算できた

