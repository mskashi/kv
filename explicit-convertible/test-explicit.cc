// test for explicit constructor

class hoge {
	public:

	double a;

	hoge() {
	}

	explicit hoge(const double& x) {
		a = x;
	}

	hoge& operator=(const double& x) {
		a = x;
		return *this;
	}
};

void func(hoge x) {
}

int main()
{
	// コンパイルできる
	hoge a(1.);

	// コンパイルできる
	hoge b = hoge(1.);

	// コンパイル出来ない。下のとの非対称性を考えるとコンパイル出来てほしい?
	// hoge c = 1.;

	// これは、代入演算子を定義したのでコンパイルできる。
	hoge d;
	d = 1.;

	// コンパイル出来ない。これはまあ自然か。
	// func(1.);

	// コンパイルできる
	func(hoge(1.));
}
