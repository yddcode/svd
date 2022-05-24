
NDouble **dmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a NDouble matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	NDouble **m;
	/* allocate pointers to rows */
	m = (NDouble **)malloc((size_t)((nrow + NR_END) * sizeof(NDouble*)));
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl] = (NDouble *)malloc((size_t)((nrow*ncol + NR_END) * sizeof(NDouble)));
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

NDouble *dvector(int nl, int nh)
/* allocate a NDouble vector with subscript range v[nl..nh] */
{
	NDouble *v;
	v = (NDouble *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(NDouble)));
	return v - nl + NR_END;
}


void getSVDuvt(NDouble *inv_rotation_mat, NDouble **u, NDouble **v, NDouble w[])
{
	/*
	 * 使用 ** 来声明，然后用dmatrix来申请空间
	 * 这样在函数引用矩阵时，就不用声明矩阵大小了
	 */
	 //NDouble **a;

	int i, j, k;
	NDouble t;
	NDouble t1[MN + 1], t2[MN + 1];
	/* 矩阵均以M,N中的最大值来申请空间，避免越界 */
	//a = dmatrix(1, MN, 1, MN);  // inv_rotation_mat

	for (i = 1; i <= M; i++) {
		for (j = 1; j <= N; j++) {
			//u[i][j] = a[i][j];

			u[i][j] = inv_rotation_mat[i * 3 + j - 4];

			//*((int *)u + i * 3 + j) = inv_rotation_mat[i * 3 + j - 4];
		}
	}

	svdcmp(u, MN, MN, w, v);

	/* Sort the singular values in descending order */
	for (i = 1; i <= N; i++) {
		for (j = i + 1; j <= N; j++) {
			if (w[i] < w[j]) { /* 对特异值排序 */
				t = w[i];
				w[i] = w[j];
				w[j] = t;

				/* 同时也要把矩阵U,V的列位置交换 */
				/* 矩阵U */
				for (k = 1; k <= M; k++) {
					t1[k] = u[k][i];
				}
				for (k = 1; k <= M; k++) {
					u[k][i] = u[k][j];
				}
				for (k = 1; k <= M; k++) {
					u[k][j] = t1[k];
				}

				/* 矩阵V */
				for (k = 1; k <= N; k++) {
					t2[k] = v[k][i];
				}
				for (k = 1; k <= N; k++) {
					v[k][i] = v[k][j];
				}
				for (k = 1; k <= N; k++) {
					v[k][j] = t2[k];
				}
			}
		}
	}


	/* 奇异值有M个，存为W矩阵，后面会用到 */
	NDouble **W;
	W = dmatrix(1, MN, 1, MN);

	for (i = 1; i <= M; i++) {
		for (j = 1; j <= N; j++) {
			if (i == j) {
				W[i][j] = w[i];
			}
			else {
				W[i][j] = 0.0;
			}
		}
	}

	/* V为NxN矩阵 */
	/* 验证结果，即计算U*W*V’看是否等于a */
	NDouble **temp;
	NDouble sum;
	temp = dmatrix(1, MN, 1, MN);
	/* 先算 U*W */
	for (i = 1; i <= M; i++) {
		for (j = 1; j <= M; j++) {
			sum = 0;
			for (k = 1; k <= N; k++) {
				sum += u[i][k] * W[k][j];
			}
			temp[i][j] = sum;
		}
	}

	/* 再算 temp*V */
	/* 先对v进行矩阵转置，存为V */
	NDouble ** V;
	V = dmatrix(1, MN, 1, MN);
	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			V[i][j] = v[j][i];
		}
	}

	return;
}


/******************************************************************************/
void svdcmp(NDouble **a, int m, int n, NDouble *w, NDouble **v)
/*******************************************************************************
Given a matrix a[1..m][1..n], this routine computes its singular value
decomposition, A = U.W.VT.  The matrix U replaces a on output.  The diagonal
matrix of singular values W is output as a vector w[1..n].  The matrix V (not
the transpose VT) is output as v[1..n][1..n].
*******************************************************************************/
{
	int flag, i, its, j, jj, k, l, nm;
	NDouble anorm, c, f, g, h, s, scale, x, y, z, *rv1;

	rv1 = dvector(1, n);
	g = scale = anorm = 0.0; /* Householder reduction to bidiagonal form */
	for (i = 1; i <= n; i++) {
		l = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m) {
			for (k = i; k <= m; k++) 
				scale += fabs(a[k][i]);
			if (scale) {
				for (k = i; k <= m; k++) 
				{
					a[k][i] /= scale;
					s += a[k][i] * a[k][i];
				}
				f = a[i][i];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][i] = f - g;
				for (j = l; j <= n; j++) 
				{
					for (s = 0.0, k = i; k <= m; k++) 
						s += a[k][i] * a[k][j];

					f = s / h;
					for (k = i; k <= m; k++) 
						a[k][j] += f * a[k][i];
				}
				for (k = i; k <= m; k++) 
					a[k][i] *= scale;
			}
		}

		w[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m && i != n) {
			for (k = l; k <= n; k++) 
				scale += fabs(a[i][k]);

			if (scale) {
				for (k = l; k <= n; k++) {
					a[i][k] /= scale;
					s += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][l] = f - g;
				for (k = l; k <= n; k++) 
					rv1[k] = a[i][k] / h;
				for (j = l; j <= m; j++) {
					for (s = 0.0, k = l; k <= n; k++) 
						s += a[j][k] * a[i][k];
					for (k = l; k <= n; k++) 
						a[j][k] += s * rv1[k];
				}
				for (k = l; k <= n; k++) 
					a[i][k] *= scale;
			}
		}
		anorm = DMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
	}
	for (i = n; i >= 1; i--) { /* Accumulation of right-hand transformations. */
		if (i < n) {
			if (g) {
				for (j = l; j <= n; j++) /* Double division to avoid possible underflow. */
					v[j][i] = (a[i][j] / a[i][l]) / g;
				for (j = l; j <= n; j++) {
					for (s = 0.0, k = l; k <= n; k++) 
						s += a[i][k] * v[k][j];
					for (k = l; k <= n; k++) 
						v[k][j] += s * v[k][i];
				}
			}
			for (j = l; j <= n; j++) 
				v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}
	for (i = IMIN(m, n); i >= 1; i--) { /* Accumulation of left-hand transformations. */
		l = i + 1;
		g = w[i];
		for (j = l; j <= n; j++) a[i][j] = 0.0;
		if (g) {
			g = 1.0 / g;
			for (j = l; j <= n; j++) {
				for (s = 0.0, k = l; k <= m; k++) 
					s += a[k][i] * a[k][j];
				f = (s / a[i][i])*g;
				for (k = i; k <= m; k++) 
					a[k][j] += f * a[k][i];
			}
			for (j = i; j <= m; j++) 
				a[j][i] *= g;
		}
		else for (j = i; j <= m; j++) 
			a[j][i] = 0.0;
		++a[i][i];
	}

	for (k = n; k >= 1; k--) { /* Diagonalization of the bidiagonal form. */
		for (its = 1; its <= 30; its++) {
			flag = 1;
			for (l = k; l >= 1; l--) { /* Test for splitting. */
				nm = l - 1; /* Note that rv1[1] is always zero. */
				if ((NDouble)(fabs(rv1[l]) + anorm) == anorm) {
					flag = 0;
					break;
				}
				if ((NDouble)(fabs(w[nm]) + anorm) == anorm) break;
			}
			if (flag) {
				c = 0.0; /* Cancellation of rv1[l], if l > 1. */
				s = 1.0;
				for (i = l; i <= k; i++) {
					f = s * rv1[i];
					rv1[i] = c * rv1[i];
					if ((NDouble)(fabs(f) + anorm) == anorm) break;
					g = w[i];
					h = pythag(f, g);
					w[i] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					for (j = 1; j <= m; j++) {
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = y * c + z * s;
						a[j][i] = z * c - y * s;
					}
				}
			}
			z = w[k];
			if (l == k) { /* Convergence. */
				if (z < 0.0) { /* Singular value is made nonnegative. */
					w[k] = -z;
					for (j = 1; j <= n; j++) 
						v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) 
				printf("no convergence in 30 svdcmp iterations\n");
			x = w[l]; /* Shift from bottom 2-by-2 minor. */
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z)*(y + z) + (g - h)*(g + h)) / (2.0*h*y);
			g = pythag(f, 1.0);
			f = ((x - z)*(x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
			c = s = 1.0; /* Next QR transformation: */
			for (j = l; j <= nm; j++) {
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;
				for (jj = 1; jj <= n; jj++) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x * c + z * s;
					v[jj][i] = z * c - x * s;
				}
				z = pythag(f, h);
				w[j] = z; /* Rotation can be arbitrary if z = 0. */
				if (z) {
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				for (jj = 1; jj <= m; jj++) {
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y * c + z * s;
					a[jj][i] = z * c - y * s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
	free_dvector(rv1, 1, n);
}


NDouble pythag(NDouble a, NDouble b)
/* compute (a2 + b2)^1/2 without destructive underflow or overflow */
{
	NDouble absa, absb;
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb) return absa * sqrt(1.0 + (absb / absa)*(absb / absa));
	else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + (absa / absb)*(absa / absb)));
}

void free_dvector(NDouble *v, int nl, int nh)
/* free a NDouble vector allocated with dvector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

#-------------------------------------
std::string file_json = file_path + "/calib_para1.json";
//std::string file_json = "D:\\Admin-1\\calib_para.json";

//将数据文件数据存入数组
	//char message[1024];

	//FILE *fp;
	//fp = fopen(file_json.c_str(), "r");
	//
	//fread(message, 1024, 1, fp);
	//cout << "\n message:: " << message << endl;

	//cJSON* cjson_test = NULL;
	//cJSON* cjson_skill = NULL;
	//int    skill_array_size = 0, i = 0;
	//cJSON* cjson_skill_item = NULL;

	///* 解析整段JSO数据 */
	//cjson_test = cJSON_Parse(message);
	//if (cjson_test == NULL) printf("parse fail.\n"); 

	/* 解析数组 */
	//cjson_skill = cJSON_GetObjectItem(cjson_test, "cameraMatrixL");
	//skill_array_size = cJSON_GetArraySize(cjson_skill);
	//printf(" skill_array_size %d \n", skill_array_size);
	//for (i = 0; i < skill_array_size; i++)
	//{
	//	cjson_skill_item = cJSON_GetArrayItem(cjson_skill, i);
	//	printf("%lf ", cjson_skill_item->valuedouble);
	//	cameraMatrixL[i] = cjson_skill_item->valuedouble;
	//}

	///* 解析数组 */
	//cjson_skill = cJSON_GetObjectItem(cjson_test, "distCoeffL");
	//skill_array_size = cJSON_GetArraySize(cjson_skill);
	//printf(" skill_array_size %d \n", skill_array_size);
	//for (i = 0; i < skill_array_size; i++)
	//{
	//	cjson_skill_item = cJSON_GetArrayItem(cjson_skill, i);
	//	printf("%lf ", cjson_skill_item->valuedouble);
	//	distCoeffL[i] = cjson_skill_item->valuedouble;
	//}

	///* 解析数组 */
	//cjson_skill = cJSON_GetObjectItem(cjson_test, "cameraMatrixR");
	//skill_array_size = cJSON_GetArraySize(cjson_skill);
	//printf(" skill_array_size %d \n", skill_array_size);
	//for (i = 0; i < skill_array_size; i++)
	//{
	//	cjson_skill_item = cJSON_GetArrayItem(cjson_skill, i);
	//	printf("%lf ", cjson_skill_item->valuedouble);
	//	cameraMatrixR[i] = cjson_skill_item->valuedouble;
	//}

	///* 解析数组 */
	//cjson_skill = cJSON_GetObjectItem(cjson_test, "distCoeffR");
	//skill_array_size = cJSON_GetArraySize(cjson_skill);
	//printf(" skill_array_size %d \n", skill_array_size);
	//for (i = 0; i < skill_array_size; i++)
	//{
	//	cjson_skill_item = cJSON_GetArrayItem(cjson_skill, i);
	//	printf("%lf ", cjson_skill_item->valuedouble);
	//	distCoeffR[i] = cjson_skill_item->valuedouble;
	//}

	///* 解析数组 */
	//cjson_skill = cJSON_GetObjectItem(cjson_test, "Rotation_mat");
	//skill_array_size = cJSON_GetArraySize(cjson_skill);
	//printf(" skill_array_size %d \n", skill_array_size);
	//for (i = 0; i < skill_array_size; i++)
	//{
	//	cjson_skill_item = cJSON_GetArrayItem(cjson_skill, i);
	//	printf("%lf ", cjson_skill_item->valuedouble);
	//	Rotation_mat[i] = cjson_skill_item->valuedouble;
	//}

	///* 解析数组 */
	//cjson_skill = cJSON_GetObjectItem(cjson_test, "Translation_vector");
	//skill_array_size = cJSON_GetArraySize(cjson_skill);
	//printf(" skill_array_size %d \n", skill_array_size);
	//for (i = 0; i < skill_array_size; i++)
	//{
	//	cjson_skill_item = cJSON_GetArrayItem(cjson_skill, i);
	//	printf("%lf ", cjson_skill_item->valuedouble);
	//	Translation_vector[i] = cjson_skill_item->valuedouble;
	//}
	
	/*cv::FileStorage fs_read(file_intrisics, cv::FileStorage::READ);
	if (fs_read.isOpened()) {
		
		fs_read["cameraMatrixL"] >> cameraMatrixL;
		std::cout << "cameraMatrixL=:" << cameraMatrixL << std::endl;
		fs_read["cameraDistcoeffL"] >> distCoeffL;
		fs_read["cameraMatrixR"] >> cameraMatrixR;
		fs_read["cameraDistcoeffR"] >> distCoeffR;
		fs_read.release();
	}
	fs_read.open(file_extrinsics, cv::FileStorage::READ);
	if (fs_read.isOpened()) {
		
		fs_read["R"] >> Rotation_mat;
		fs_read["T"] >> Translation_vector;
		fs_read.release();
	}*/
  
  
  //recPara = { { 4.641904159148351e+02, 0.095917024363799, 3.296733848523218e+02, 0, 4.639347248969304e+02, 2.467155339880595e+02, 0, 0, 1},
	//		{-0.005740358815623, 0.150138928440197, 2.502755501700755e-04, -1.534733915894606e-05, -0.365559238123493},
	//		{4.641168066978353e+02, -0.174655912482272, 3.120078587449918e+02, 0, 4.639998738833286e+02, 2.342449563713411e+02, 0, 0, 1},
	//		{-0.010717391774093, 0.125566848760184, -3.888266645260528e-04, 3.162241263086794e-04, -0.311054863206094},
	//		640, 480, {0.999949550593789, 3.621319138589388e-04, -0.010038183488850, -2.800018093138742e-04, 0.999966487938929, 0.008181967860471, 0.010040810040312, -0.008178744375512, 0.999916141620973},
	//		{-17.486448195277198, -0.103882128996959, -0.735345403353312} };

	//NDouble cameraMatrixL[3 * 3] = { 4.641904159148351e+02, 0.095917024363799, 3.296733848523218e+02, 0, 4.639347248969304e+02, 2.467155339880595e+02, 0, 0, 1 };
	//NDouble distCoeffL[5] = { -0.005740358815623, 0.150138928440197, 2.502755501700755e-04, -1.534733915894606e-05, -0.365559238123493 };
	//NDouble cameraMatrixR[3 * 3] = { 4.641168066978353e+02, -0.174655912482272, 3.120078587449918e+02, 0, 4.639998738833286e+02, 2.342449563713411e+02, 0, 0, 1 };
	//NDouble distCoeffR[5] = { -0.010717391774093, 0.125566848760184, -3.888266645260528e-04, 3.162241263086794e-04, -0.311054863206094 };
	//NDouble Rotation_mat[3 * 3] = { 0.999949550593789, 3.621319138589388e-04, -0.010038183488850, 
	//								-2.800018093138742e-04, 0.999966487938929, 0.008181967860471, 
	//								0.010040810040312, -0.008178744375512, 0.999916141620973 };
	///*NDouble Rotation_mat[3 * 3] = { 0.999949550593789, -2.800018093138742e-04, 0.010040810040312,
	//								3.621319138589388e-04, 0.999966487938929, -0.008178744375512,
	//								-0.010038183488850, 0.008181967860471, 0.999916141620973 };*/
	//NDouble Translation_vector[3] = { -17.486448195277198, -0.103882128996959, -0.735345403353312 };
	
	// camera 2	
	//NDouble cameraMatrixL[3 * 3] = { 4.737129766937609e+02, -0.426031352804459, 3.460542641062556e+02, 0, 4.734087347513511e+02, 2.321633553731306e+02, 0, 0, 1 };
	//NDouble distCoeffL[5] = { 6.025743531601421e-04, 0.061064553149077, -0.001580623325185, -8.434739176894420e-04, -0.159862909985620 };
	//NDouble cameraMatrixR[3 * 3] = { 4.708238203143260e+02, -0.395343454006348, 3.401949131777267e+02, 0, 4.704118217638651e+02, 2.333316212724313e+02, 0, 0, 1 };
	//NDouble distCoeffR[5] = { -0.009935311490590, 0.064838823786173, -0.001022943456134, -2.254687497919127e-04, -0.141943102227453 };
	//NDouble Rotation_mat[3 * 3] = { 0.999999786801833, 3.689293258164731e-05, -6.519472373428515e-04, 
	//								-3.929587731533848e-05, 0.999993205287634, -0.003686167982794, 
	//								6.518068140020233e-04, 0.003686192815749, 0.999992993540656 };
	//NDouble Translation_vector[3] = { -17.019522869698040, -0.058710480933973, 0.367647336599477 };

	// camera RGB 
	/*NDouble cameraMatrixL[3 * 3] = { 1.326671630149860e+03, -0.016745645393437, 7.669422054166025e+02, 0, 1.325731837380084e+03, 7.223675450206308e+02, 0, 0, 1 };
	NDouble distCoeffL[5] = { 0.047403550290918, -0.045262497163565, -2.177909471803504e-04, -3.239241591594921e-05, 0.010497557606732 };
	NDouble cameraMatrixR[3 * 3] = { 1.329276298464530e+03, 0.347520822627361, 7.843901998286794e+02, 0, 1.328604588922002e+03, 7.725371097358284e+02, 0, 0, 1 };
	NDouble distCoeffR[5] = { 0.046561647354570, -0.033461777040866, 5.020081333798824e-04, 8.348591288957536e-04, -0.015523326186209 };
	NDouble Rotation_mat[3 * 3] = { 0.999952809527407,  -4.599352010183297e-05, -0.009714762109355,
									8.763264419291058e-05, 0.999990812128495, 0.004285788038768,
									0.009714475732891, -0.004286437120697, 0.999943626120016 };
	NDouble Translation_vector[3] = { -66.771585258455810, 0.067907235516853, 0.366560509710477 };*/

	// camera RGB 1520
	/*NDouble cameraMatrixL[3 * 3] = { 1.325284349267375e+03, 0, 7.676686952332727e+02, 0, 1.324292528888089e+03, 7.217757152748692e+02, 0, 0, 1 };
	NDouble distCoeffL[5] = { 0.045834970552427, -0.043067850920675, -4.864790422435129e-04, 1.100836601500052e-04, 0.0103230393596036 };
	NDouble cameraMatrixR[3 * 3] = { 1.335618983983285e+03, 0, 7.813730292987332e+02, 0, 1.335208737157443e+03, 7.727521086311323e+02, 0, 0, 1 };
	NDouble distCoeffR[5] = { 0.044656995422689, -0.028677035147351, 2.285277881159679e-04, 6.907043137355630e-04, -0.020893293064865 };
	NDouble Rotation_mat[3 * 3] = { 0.999947052449964, -1.266436996468097e-04, -0.010289618943473,
									1.709601365275807e-04, 0.999989771082182, 0.004306142066887,
									0.010288978047418, -0.004307673181874, 0.999950074706269 };
	NDouble Translation_vector[3] = { -66.821333030210200, 0.04091818266155, 0.30573833842317 };*/

	// camera RGB 1520 2
	/*NDouble cameraMatrixL[3 * 3] = { 1.3275677813985833e+03, 0., 7.6755327715067710e+02, 0., 1.3265838745708279e+03, 7.2101758587381835e+02, 0., 0., 1. };
	NDouble distCoeffL[5] = { 4.8594796161606203e-02, -3.2411265932953959e-02,	-3.0068218777243224e-04, 2.4708039409647572e-04,	-1.3420603209160562e-02 };
	NDouble cameraMatrixR[3 * 3] = { 1.3286547807928778e+03, 0., 7.8233618878204959e+02, 0.,	1.3276137188018470e+03, 7.7155109021414546e+02, 0., 0., 1. };
	NDouble distCoeffR[5] = { 5.1133114529946781e-02, -4.6139155743770310e-02,	3.7140059965429896e-04, 5.2987953504707991e-04,	-5.0557456485294016e-03 };
	NDouble Rotation_mat[3 * 3] = { 9.9993340197198310e-01, -1.3581936124299688e-05,	1.1540859424993839e-02, 6.7734408651969481e-05,
	9.9998899086445125e-01, -4.6918612454242524e-03,	-1.1540668645548341e-02, 4.6923304900060650e-03,	9.9992239448958553e-01 };
	NDouble Translation_vector[3] = { -6.6668083361254361e+01, 1.9031883916356995e-01,	-3.6524448819466798e-01 };*/
	// 齐感
	/*NDouble cameraMatrixL[3 * 3] = { 3.795368131816281e+02, 0, 3.237930819565501e+02, 0, 3.795370175656677e+02, 2.470746434790044e+02, 0, 0, 1 };
	NDouble distCoeffL[5] = { -0.123920866431494, 0.027265191142467, 0.002108190741892, 3.026942273651654e-04, 0.335307211332364 };
	NDouble cameraMatrixR[3 * 3] = { 3.792158891776506e+02, 0, 3.105319740260449e+02, 0, 3.791527957714079e+02, 2.576409190398907e+02, 0, 0, 1 };
	NDouble distCoeffR[5] = { -0.117425915591507, -0.038835862377629, 4.643801957269621e-04, -7.913619172457312e-04, 0.566395769423928 };
	NDouble Rotation_mat[3 * 3] = { 9.9987522622424430e-01, 1.548905465397049e-04,-0.001902054144030,
			 -1.675114503289458e-04,9.9999679087215776e-01, -0.006637041593578,	
			0.001900984208835, 0.006637348124023,9.9987305678811378e-01 };
	NDouble Translation_vector[3] = { -0.306227509126052, -17.616441563936630, 2.010404007655490 };*/
	// qigan shuiping xuanzhuan
	/*NDouble cameraMatrixL[3 * 3] = { 3.787963649257394e+02, 0, 2.447480635608792e+02, 0, 3.788802314769000e+02, 3.206864381507143e+02, 0, 0, 1 };
	NDouble distCoeffL[5] = { -0.128915672023569, 0.149908196169193, 0.001677882134990, 2.409266583331020e-04, 0.052996542741338 };
	NDouble cameraMatrixR[3 * 3] = { 3.775344249731961e+02, 0, 2.597564967613446e+02, 0, 3.775720373393367e+02, 3.330289396957376e+02, 0, 0, 1 };
	NDouble distCoeffR[5] = { -0.134205374998256, 0.165914586173352, 0.002669413474095, 0.001862156332013, -0.077410879103077 };
	NDouble Rotation_mat[3 * 3] = { 9.9987522622424430e-01, 0.001673056586856,0.004420038918792,
			 -0.001671831396112,9.9999679087215776e-01, -2.808702958918097e-04,
			-0.004420502479293, 2.734775992976784e-04,9.9987305678811378e-01 };
	NDouble Translation_vector[3] = { -17.653423775793424, 0.155062642959844, 0.393790941494345 };*/
	// qigan opencv +c+	

	/*NDouble cameraMatrixL[3 * 3] = { 3.8113711696578571e+02, 0., 2.4619028318811866e+02, 0.,
	   3.8119176049104567e+02, 3.1803467056198241e+02, 0., 0., 1. };
	NDouble distCoeffL[5] = { -1.4099573926588527e-01, 2.0599252528048598e-01,
	   1.0593704242486961e-03, 1.4673690631199548e-03,
	   -1.2686559843569969e-01 };
	NDouble cameraMatrixR[3 * 3] = { 3.8013460965888106e+02, 0., 2.5750734582361684e+02, 0.,
	   3.7984241880450639e+02, 3.3199685002428799e+02, 0., 0., 1. };
	NDouble distCoeffR[5] = { -1.3657634501952201e-01, 1.5664704319424846e-01,
	   3.4474936385531954e-03, 3.3537232750499064e-04,
	   -8.6463781433560885e-02 };
	NDouble Rotation_mat[3 * 3] = { 9.9998591253534874e-01, -1.3039180236689371e-03,
		5.1453404778846659e-03, 1.3226479154873733e-03,
		9.9999250613773727e-01, -3.6384434665198886e-03,
		-5.1405576873977416e-03, 3.6451976839337059e-03,
		9.9998014340311159e-01 };
	NDouble Translation_vector[3] = { -1.7470836313886846e+01, -8.0145458423687685e-02,	   8.3857200972890500e-01 };*/
  
  
