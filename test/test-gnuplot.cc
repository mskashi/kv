#include <kv/gnuplot.hpp>

int main()
{
	kv::gnuplot g;

	// 初期化
	g.open("/usr/bin/gnuplot");
	// ファイルに吐くには例えばこちら
	// g.open("/usr/bin/gnuplot", "postscript eps color", "aaa.eps");
	// teeを使って描画コマンドをファイルに残す方法もある。
	// g.open("tee aaa.dat | /usr/bin/gnuplot");

	// 描画範囲設定
	g.screen(0, 0, 10, 10);

	g.line(1,1,3,4);
	g.point(4,4);
	g.ellipsef(8,2,2,1);
	g.rect(5,5,6,7);
	g.rectf(8,8,9,9);

	getchar();

	// 描画終了
	g.close();
}
