<pre xml:lang="latex">\sqrt{2}</pre>
```math
SE = \frac{\sigma}{\sqrt{n}}
```

# 010TwoBody

分子動力学法の最も単純な例として、質量1の2つの質点がバネ定数1のバネでつながれた系の
運動(座標の時間変化)を、運動方程式を数値的に解くことで求めます。

物体の運動は、ニュートンの運動方程式

<img src="https://latex.codecogs.com/svg.image?\begin{equation}  {\bf F} = m {\bf a}\end{equation}" />

に従います。ここで、Fはバネが物体に及ぼす力、mは質量、aは加速度です。

バネの力は、2つの物体の相対座標で決まります。また、加速度は座標の時間での2回微分ですから、この式を物体のかたほうについてもっと具体的に書くと、

<img src="https://latex.codecogs.com/svg.image?\begin{equation}{\bf F}_1 = -k ({\bf r}_1- {\bf r}_2) = m_1 {\bf a}_1 \end{equation}" />

と書けます。力$\bf F$と位置$\bf r$と加速度$\bf a$はベクトルです。$k$はバネ定数です．作用反作用の法則により，もう一方の末端には逆向きの力が加わります．

These equations can also be derived from the potential energy function of a spring by differentiation.

<img src="https://latex.codecogs.com/svg.image?
\begin{matrix}
V({\bf r}_1, {\bf r}_2) &=& \frac{k}{2}\left({\bf r}_1-{\bf    r}_2\right)^2\\
 {\bf F}_1&=& -\frac{\partial V}{\partial{\bf r}_1}\\
   &=&-k\left({\bf r}_1-{\bf r}_2\right)\\
{\bf F}_2&=& -\frac{\partial V}{\partial{\bf r}_2}\\
   &=&-{\bf F}_1 \end{matrix}
" />

ところで、加速度は速度を時間で微分したもの

<img src="https://latex.codecogs.com/svg.image?\begin{equation}a_x = \partial v_x / \partial t\end{equation}" />

また速度は座標を時間で微分したものです。

<img src="https://latex.codecogs.com/svg.image?\begin{equation}v_x = \partial r_x / \partial t\end{equation}" />

$_x$は速度や加速度ベクトルの$x$成分を表します．$y$，$z$成分も同様です．

これを、差分で近似的に表すと、(添字$_x$は省略しました)

<img src="https://latex.codecogs.com/svg.image?\begin{matrix}a(t)&=&[ v(t+\Delta t) - v(t) ] / \Delta t\\v(t)&=&[ r(t+\Delta t) - r(t) ] / \Delta t\end{matrix}" />

移項して書きなおすと、

<img src="https://latex.codecogs.com/svg.image?
\begin{matrix}
v(t+\Delta t) = v(t) + a(t) \Delta t\\
r(t+\Delta t) = r(t) + v(t) \Delta t
\end{eqnarray}" />

と書けます。つまり、現在の座標から力を計算し、力から運動方程式によって現在の加速度を計算し、現在の速度と加速度から少し未来の速度を、そして現在の座標と速度から少し未来の座標を計算できることがわかります。

この計算を繰り返し行うことで、物体の未来の位置と速度を計算できます。これが分子動力学法の中核です。プログラム`main.f90`を上の式と見比べて、どこでどの計算をやっているかを想像してみて下さい。

このプログラムを実行するには、コンピュータが理解できる形に「翻訳」する必要があります。この作業をコンパイルと呼びます。`main.f90`をコンパイルして、実行ファイル`main`を作るには、ターミナルで

`gfortran main.f90 -o main`

と入力します。

エラーなくコンパイルできたら、

`./main`

と入力して下さい。`./`は、"現在のディレクトリにある"、という意味を表しています。画面には何か数字がたくさん表示され、いかにもたくさん計算している雰囲気になります。

ちゃんとプログラムが動いていて、物理法則にのっとった運動をしているかどうかを確認するには、エネルギー保存則をチェックするのが簡単です。この系の場合、ポテンシャルエネルギーはバネのエネルギー

<img src="https://latex.codecogs.com/svg.latex?
\begin{equation}\label{eq:e_sp}
E_p = k \left| r_1 - r_2 \right|^2 / 2
\end{equation}
" />

運動エネルギーは2つの物体それぞれ

<img src="https://latex.codecogs.com/svg.image?
\begin{equation}
E_k = m \left| v \right|^2 / 2
\end{equation}
" />

です。プログラムでは、これらの値を画面に出力するようになっています。

数値計算ですから、全エネルギーが厳密に一定になるわけではありませんが、運動エネルギーやポテンシャルエネルギーの時間変動に比べると、それらの和はほとんど変動しないことを、`gnuplot`などでプロットして確認して下さい。

## 練習問題1

バネのエネルギーを末端間距離$r$で微分すると，バネの両端に加わる力が求められます．

<img src="https://latex.codecogs.com/svg.image?
\begin{matrix}
E_p(r)&=&kr^2/2\\
{\bf F}&=&-\frac{\partial E_p(r)}{\partial {\bf r}}\\
&=&-\frac{\partial E_p(r)}{\partial r}\frac{\partial r}{\partial {\bf r}}\\
&=&-k{\bf r}
\end{matrix}
" />

クーロンポテンシャル関数$V(r) = k / r$を微分して、クーロンポテンシャルで相互作用する2つの物体に加わる力のベクトルを数式として求めて下さい。


## 練習問題2

実際に2つの物体(2次元でも3次元でも構わない。座標はてきとうに決める)の間に働くポテンシャルエネルギーと、物体の座標をわずかにずらした時のエネルギー変化から、上の式が正しいかどうかを検証するプログラムを書きましょう。

