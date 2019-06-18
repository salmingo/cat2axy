/*
 Name        : cat2axy.cpp
 Author      : Xiaomeng Lu
 Copyright   : SVOM@NAOC, CAS
 Description : 将SExtractor生成的cat文件转换为astrometry.net的axy文件
 目标: 将sample.cat转换为sample.axy
 (1) 按照流量对cat文件数据重新排序. cat文件各列对应:
 Col 1: X
 Col 2: Y
 Col 3: Flux
 Col 4: FWHM
 Col 5: Elongation
 (2) 将图像划分为小区块, 每个区块最多包含3颗星
 (3) 将星象X/Y坐标写入axy文件

 @version 0.1
 @date 2019-06-08
 @note
 - 约束条件: 恒星延展率小于2.0
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <longnam.h>
#include <fitsio.h>
#include <boost/filesystem.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/make_shared.hpp>

using namespace std;

#define GRID		128
#define NDX_X		0
#define NDX_Y		1
#define NDX_FLUX	2
#define NDX_FWHM	3
#define NDX_ELON	4
#define NDX_NUM		(NDX_ELON + 1)

struct aobject {
	float info[NDX_NUM];
};
typedef boost::shared_ptr<aobject> objptr;
typedef vector<objptr> objptrvec;

bool less_brightness(objptr x1, objptr x2) {
	return (x1->info[NDX_FLUX] > x2->info[NDX_FLUX]);
}

/*!
 * @brief 从cat文件中加载星象信息
 * @param filepath
 * @param objs
 */
bool load_cat(const char *filepath, objptrvec &objs) {
	// 加载数据
	int pos;
	char line[200], seps[] = " \r", *token;
	float snr;
	FILE *fpcat = fopen(filepath, "r");

	if (fpcat) {
		while (!feof(fpcat)) {
			if (!fgets(line, 200, fpcat) || line[0] == '#')
				continue;

			pos = -1;
			token = strtok(line, seps);
			objptr one = boost::make_shared<aobject>();
			while (token && ++pos < NDX_NUM) {
				one->info[pos] = atof(token);
				token = strtok(NULL, seps);
			}
			if (one->info[NDX_FLUX] > 30.0 && one->info[NDX_FWHM] > 1.0
					&& one->info[NDX_ELON] < 2.0)
				objs.push_back(one);
		}
		fclose(fpcat);

		// 按照流量排序
		sort(objs.begin(), objs.end(), less_brightness);
		return true;
	}
	return false;
}

void select_refstar(objptrvec& objs, int w, int h, objptrvec& refs) {
	int nw = (w + GRID - 1) / GRID;
	int nh = (h + GRID - 1) / GRID;

	if (nw * nh < 4) {
		for (objptrvec::iterator it = objs.begin(); it != objs.end(); ++it)
			refs.push_back(*it);
	}
	else {
		int x0 = (w % GRID) / 2;
		int y0 = (h % GRID) / 2;
		int i, j, k;
		char *map = new char[nw * nh];
		memset(map, 0, nw * nh);

		for (objptrvec::iterator it = objs.begin(); it != objs.end(); ++it) {
			i = ((*it)->info[NDX_X] - x0) / GRID;
			j = ((*it)->info[NDX_Y] - y0) / GRID;
			if (i == nw || j == nh)
				continue;
			if (map[k = i * j] <= 5) {
				++map[k];
				refs.push_back(*it);
			}
		}
		delete[] map;
	}
}

void output_axy(objptrvec& refs, const char *filepath) {
	namespace fs = boost::filesystem;
	fs::path path = filepath;
	if (fs::exists(path))
		fs::remove(path);
	// 创建文件并写入数据
	fitsfile *fitsptr;	//< 基于cfitsio接口的文件操作接口
	int status(0);
	char *ttype[10] = { "X", "Y" };
	char *tform[10] = { "E", "E" };
	int n = refs.size();
	float *x = new float[n];
	float *y = new float[n];

	for (int i = 0; i < n; ++i) {
		x[i] = refs[i]->info[NDX_X];
		y[i] = refs[i]->info[NDX_Y];
	}

	fits_create_file(&fitsptr, filepath, &status);
	fits_insert_btbl(fitsptr, refs.size(), 2, ttype, tform, NULL, NULL, 0,
			&status);
	fits_write_col(fitsptr, TFLOAT, 1, 1, 1, n, x, &status);
	fits_write_col(fitsptr, TFLOAT, 2, 1, 1, n, y, &status);
	fits_close_file(fitsptr, &status);
	delete[] x;
	delete[] y;

	if (status) {
		char txt[100];
		ffgerr(status, txt);
		cout << txt << endl;
	}
}

/*
 * 命令行参数:
 * 1 - cat文件路径
 * 2 - 图像宽度
 * 3 - 图像高度
 */
int main(int argc, char **argv) {
	if (argc < 4) {
		cout << "Usage:" << endl;
		cout << "\t cat2axy CAT_Path X Y" << endl;
		for (int i = 0; i < argc; ++i)
			cout << argv[i] << "  ";
		cout << endl;
		return -1;
	}

	namespace fs = boost::filesystem;
	fs::path pathcat = argv[1];
	fs::path pathaxy = pathcat;
	int w = atoi(argv[2]);
	int h = atoi(argv[3]);
	objptrvec objs, refs;

	pathaxy.replace_extension(fs::path(".axy"));
	if (load_cat(pathcat.c_str(), objs)) {
		select_refstar(objs, w, h, refs);
		if (refs.size() >= 5)
			output_axy(refs, pathaxy.c_str());
		else
			cout << "not found enough reference stars" << endl;
		objs.clear();
		refs.clear();
	}
	else
		cout << "failed to load CAT file: " << argv[1] << endl;

	return 0;
}
