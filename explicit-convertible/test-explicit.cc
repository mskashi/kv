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
	// ����ѥ���Ǥ���
	hoge a(1.);

	// ����ѥ���Ǥ���
	hoge b = hoge(1.);

	// ����ѥ������ʤ������ΤȤ����о�����ͤ���ȥ���ѥ������Ƥۤ���?
	// hoge c = 1.;

	// ����ϡ������黻�Ҥ���������Τǥ���ѥ���Ǥ��롣
	hoge d;
	d = 1.;

	// ����ѥ������ʤ�������Ϥޤ���������
	// func(1.);

	// ����ѥ���Ǥ���
	func(hoge(1.));
}
