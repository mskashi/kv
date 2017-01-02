/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef GNUPLOT_HPP
#define GNUPLOT_HPP

#include <stdio.h>

#ifdef _WIN32
#define popen _popen
#define pclose _pclose
#define DEFAULT_GNUPLOT "\"C:\\Program Files (x86)\\gnuplot\\bin\\gnuplot.exe\""
#else
#define DEFAULT_GNUPLOT "gnuplot"
#endif

namespace kv {

class gnuplot {
	public:

	FILE *p;

	bool open(const char *s=DEFAULT_GNUPLOT, const char *term=NULL, const char *output=NULL) {
		p = popen(s, "w");
		if (p == NULL) return false;
		if (term != NULL) {
			fprintf(p, "set terminal %s\n", term);
			fflush(p);
		}
		if (output != NULL) {
			fprintf(p, "set output \"%s\"\n", output);
			fflush(p);
		}
		send_command("set multiplot");
		return true;
	}

	bool close() {
		send_command("quit");
		if (pclose(p) == -1) return false;
		return true;
	}

	void screen(double x1, double y1, double x2, double y2) const {
		fprintf(p, "set xrange [%.17f:%.17f]\n", x1, x2);
		fprintf(p, "set yrange [%.17f:%.17f]\n", y1, y2);
		fflush(p);
	}

	void line(double x1, double y1, double x2, double y2, int t=1) const {
		fprintf(p, "plot '-' notitle with lines linetype %d\n", t);
		fprintf(p, "%.17f %.17f\n", x1, y1);
		fprintf(p, "%.17f %.17f\n", x2, y2);
		#if 0
		fprintf(p, "plot '-' notitle with vectors nohead lt %d\n", t);
		fprintf(p, "%.17f %.17f %.17f %.17f\n", x1, y1, x2-x1, y2-y1);
		#endif
		fprintf(p, "end\n");
		fflush(p);
	}

	void point(double x, double y, int t=6) const {
		fprintf(p, "plot '-' notitle with points linetype %d\n", t);
		fprintf(p, "%.17f %.17f\n", x, y);
		fprintf(p, "end\n");
		fflush(p);
	}

	void rect(double x1, double y1, double x2, double y2, int t=1) const {
		fprintf(p, "plot '-' notitle with lines linetype %d\n", t);
		fprintf(p, "%.17f %.17f\n", x1, y1);
		fprintf(p, "%.17f %.17f\n", x1, y2);
		fprintf(p, "%.17f %.17f\n", x2, y2);
		fprintf(p, "%.17f %.17f\n", x2, y1);
		fprintf(p, "%.17f %.17f\n", x1, y1);
		fprintf(p, "end\n");
		fflush(p);
	}

	void rectf(double x1, double y1, double x2, double y2, int t=2) const {
		fprintf(p, "plot '-' notitle with filledcurves linetype %d\n", t);
		fprintf(p, "%.17f %.17f\n", x1, y1);
		fprintf(p, "%.17f %.17f\n", x1, y2);
		fprintf(p, "%.17f %.17f\n", x2, y2);
		fprintf(p, "%.17f %.17f\n", x2, y1);
		fprintf(p, "%.17f %.17f\n", x1, y1);
		fprintf(p, "end\n");
		fflush(p);
	}

	void ellipse(double cx, double cy, double rx, double ry, int t=1) const {
		fprintf(p, "set parametric\n");
		fprintf(p, "plot [0:2*pi] %.17f+%.17f*cos(t), %.17f+%.17f*sin(t) notitle with lines linetype %d\n", cx, rx, cy, ry, t);
		fprintf(p, "unset parametric\n");
		fflush(p);
	}

	void circle(double cx, double cy, double r, double ry, int t=1) const {
		ellipse(cx, cy, r, r, t);
	}

	void ellipsef(double cx, double cy, double rx, double ry, int t=2) const {
		fprintf(p, "set parametric\n");
		fprintf(p, "plot [0:2*pi] %.17f+%.17f*cos(t), %.17f+%.17f*sin(t) notitle with filledcurves linetype %d\n", cx, rx, cy, ry, t);
		fprintf(p, "unset parametric\n");
		fflush(p);
	}

	void circlef(double cx, double cy, double r, double ry, int t=2) const {
		ellipsef(cx, cy, r, r, t);
	}

	void send_command(const char *s) const {
		fprintf(p, "%s\n", s);
		fflush(p);
	}

	void clear() const {
		send_command("clear");
	}
};

} // namespace kv

#endif // GNUPLOT_HPP
