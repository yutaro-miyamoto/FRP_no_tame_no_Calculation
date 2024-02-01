// carbon.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include <stdio.h>
#include<math.h>
#define M_PI        3.14159265358979323846264338327950288 
#define N 3 /* 逆行列を求める行列の行数・列数 */
#define MAX_ERR 1e-10 /* 許容する誤差 */


void Display(double Q[3][3]);
int check(double mat[N][N], double inv[N][N]);
double Syokika(double Q[3][3]);
double Tenchi(double Q[3][3]);
double GyakuGyoretsu(double Q[3][3]);
double Gyoretsu_Kakesan(double Q[3][3], double R[3][3], double A[3][3]);
double Input(double A[3][3]); 
double Stress_Strain(double T_sigma[3][3], double Q[3][3], double T_epsilon[3][3], double result[3][3]);
double In_plane_stiffness_matrix_function(double Q_bar[3][3], double Thickness, double In_plane_stiffness_matrix[3][3]);
double Tensile_bending_coupling_matrix_function(double Q_bar[3][3], double Thickness, double Tensile_bending_coupling_matrix[3][3]);
double Out_of_plane_stiffness_matrix_function(double Q_bar[3][3], double Thickness, double Out_of_plane_stiffness_matrix[3][3]);

int main()
{
	int n, i;
	double E_L, E_T, nu_LT, nu_TL, G_LT, Q[3][3], T_sigma[3][3], T_epsilon[3][3], theta, h[100];
	double In_plane_stiffness_matrix[3][3], Tensile_bending_coupling_matrix[3][3], Out_of_plane_stiffness_matrix[3][3];
	Syokika(In_plane_stiffness_matrix);
	Syokika(Tensile_bending_coupling_matrix);
	Syokika(Out_of_plane_stiffness_matrix);

	printf("何層の積層を扱うか？：");
	scanf_s("%d", &n);

	for (i = 0; i < n; i++)
	{

		printf("Theta: ");
		scanf_s("%lf", &theta);
		theta = theta * M_PI / 180;
		/*θだけ回転した座標系への応力の変換行列*/
		T_sigma[0][0] = pow(cos(theta), 2);
		T_sigma[0][1] = pow(sin(theta), 2);
		T_sigma[0][2] = sin(2 * theta);
		T_sigma[1][0] = pow(sin(theta), 2);
		T_sigma[1][1] = pow(cos(theta), 2);
		T_sigma[1][2] = -sin(2 * theta);
		T_sigma[2][0] = -0.5 * sin(2 * theta);
		T_sigma[2][1] = 0.5 * sin(2 * theta);
		T_sigma[2][2] = cos(2 * theta);

		/*θだけ回転した座標系へのひずみの変換行列*/
		T_epsilon[0][0] = pow(cos(theta), 2);
		T_epsilon[0][1] = pow(sin(theta), 2);
		T_epsilon[0][2] = 0.5 * sin(2 * theta);
		T_epsilon[1][0] = pow(sin(theta), 2);
		T_epsilon[1][1] = pow(cos(theta), 2);
		T_epsilon[1][2] = -0.5 * sin(2 * theta);
		T_epsilon[2][0] = -sin(2 * theta);
		T_epsilon[2][1] = sin(2 * theta);
		T_epsilon[2][2] = cos(2 * theta);

		/*剛性行列[Q]を生成*/
		Syokika(Q);
		printf("E_L(繊維方向の弾性係数): ");//繊維方向の弾性係数
		scanf_s("%lf", &E_L);
		printf("E_T(繊維直行方向の弾性係数): ");//繊維直行方向の弾性係数
		scanf_s("%lf", &E_T);
		printf("nu_LT(主ポアソン比): ");//主ポアソン比
		scanf_s("%lf", &nu_LT);
		printf("nu_TL(従ポアソン比): ");//従ポアソン比
		scanf_s("%lf", &nu_TL);
		printf("G_LT(せん断弾性率): ");//せん断弾性率
		scanf_s("%lf", &G_LT);
		Q[0][0] = E_L / (1 - (nu_LT * nu_TL));
		Q[0][1] = nu_LT * E_T / (1 - (nu_LT * nu_TL));
		Q[1][0] = nu_TL * E_L / (1 - (nu_LT * nu_TL));
		Q[1][1] = E_T / (1 - (nu_LT * nu_TL));
		Q[2][2] = G_LT;
		/*x-y座標系の応力とひずみの関係式の係数行列Q_barを求める*/
		double Q_bar[3][3];
		Syokika(Q_bar);
		Stress_Strain(T_sigma, Q, T_epsilon, Q_bar);



		printf("%d層のh(積層板の厚さ): ", i + 1);//積層板の厚さ
		scanf_s("%lf", &h[i]);
		/*面内剛性行列(In_plane_stiffness_matrix)を求める*/
		In_plane_stiffness_matrix_function(Q_bar, h[i], In_plane_stiffness_matrix);

		/*引っ張り-曲げカップリング行列を求める*/
		Tensile_bending_coupling_matrix_function(Q_bar, h[i], Tensile_bending_coupling_matrix);

		/*面外剛性行列*/
		Out_of_plane_stiffness_matrix_function(Q_bar, h[i], Out_of_plane_stiffness_matrix);

	}


	printf("面内剛性行列[A]を表示\n=====================\n");
	Display(In_plane_stiffness_matrix);
	printf("=====================\n");
	printf("引っ張り-曲げカップリング行列[B]を表示\n=====================\n");
	Display(Tensile_bending_coupling_matrix);
	printf("=====================\n");
	printf("面外剛性行列[D]を表示\n=====================\n");
	Display(Out_of_plane_stiffness_matrix);
	printf("=====================\n");
}

void Display(double Q[3][3])//行列の表示
{
	// 表示
	int i, j;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			printf("%lf ", Q[i][j]);
		}
		printf("\n");
	}
}

int check(double mat[N][N], double inv[N][N]) {

	double inner_product;
	int i, j, k;
	double ans;
	double diff;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			inner_product = 0;
			for (k = 0; k < N; k++) {
				inner_product += mat[i][k] * inv[k][j];
			}

			ans = (i == j) ? 1 : 0;
			diff = ans - inner_product;
			if (fabs(diff) > MAX_ERR) {
				return 0;
			}
		}
	}

	return 1;
}

double Tenchi(double Q[3][3])
{
	int i, j;
	double Q_a[3][3];
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			Q_a[i][j] = Q[i][j];
		}
	}
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			Q[j][i] = Q_a[i][j];
		}
	}
	return Q[3][3];
}

double GyakuGyoretsu(double Q[3][3])
{
	/* 逆行列を求める行列用の２次元配列 */
	double mat[N][N];

	/* 逆行列用の２次元配列 */
	double inv[N][N];

	/* 掃き出し法を行う行列 */
	double sweep[N][N * 2];

	int i; /* 行 */
	int j; /* 列 */
	int k; /* 注目対角成分が存在する列 */

	double a; /* 定数倍用 */

	/* 正方行列の各成分をセット */
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			mat[i][j] = Q[i][j];
		}
	}

	for (i = 0; i < N; i++) 
	{
		for (j = 0; j < N; j++)
		{
			/* sweepの左側に逆行列を求める行列をセット */
			sweep[i][j] = mat[i][j];

			/* sweepの右側に単位行列をセット */
			sweep[i][N + j] = (i == j) ? 1 : 0;
		}
	}


	/* 全ての列の対角成分に対する繰り返し */
	for (k = 0; k < N; k++) 
	{

		/* 最大の絶対値を注目対角成分の絶対値と仮定 */
		double max = fabs(sweep[k][k]);
		int max_i = k;

		/* k列目が最大の絶対値となる行を探す */
		for (i = k + 1; i < N; i++) 
		{
			if (fabs(sweep[i][k]) > max) 
			{
				max = fabs(sweep[i][k]);
				max_i = i;
			}
		}

		if (fabs(sweep[max_i][k]) <= MAX_ERR) 
		{
			/* 逆行列は求められない */
			printf("逆行列は求められません...\n");
			return 0;
		}

		/* 操作（１）：k行目とmax_i行目を入れ替える */
		if (k != max_i) {
			for (j = 0; j < N * 2; j++)
			{
				double tmp = sweep[max_i][j];
				sweep[max_i][j] = sweep[k][j];
				sweep[k][j] = tmp;
			}
		}

		/* sweep[k][k]に掛けると1になる値を求める */
		a = 1 / sweep[k][k];

		/* 操作（２）：k行目をa倍する */
		for (j = 0; j < N * 2; j++)
		{
			/* これによりsweep[k][k]が1になる */
			sweep[k][j] *= a;
		}

		/* 操作（３）によりk行目以外の行のk列目を0にする */
		for (i = 0; i < N; i++) 
		{
			if (i == k) {
				/* k行目はそのまま */
				continue;
			}

			/* k行目に掛ける値を求める */
			a = -sweep[i][k];

			for (j = 0; j < N * 2; j++) 
			{
				/* i行目にk行目をa倍した行を足す */
				/* これによりsweep[i][k]が0になる */
				sweep[i][j] += sweep[k][j] * a;
			}
		}
	}

	/* sweepの右半分がmatの逆行列 */
	for (i = 0; i < N; i++) 
	{
		for (j = 0; j < N; j++) 
		{
			inv[i][j] = sweep[i][N + j];
		}
	}

	/*Qにinvを代入*/
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			Q[i][j] = inv[i][j];
		}
	}

	/* 逆行列invを表示 
	for (i = 0; i < N; i++) 
	{
		for (j = 0; j < N; j++) 
		{
			printf("%f, ",Q[i][j]);
		}
		printf("\n");
	}*/

	/* 検算 
	if (check(mat, inv)) 
	{
		printf("invはmatの逆行列です！！\n");
	}
	else 
	{
		printf("invはmatの逆行列になってません...\n");
	}*/

	return Q[3][3];
}

double Syokika(double Q[3][3])//行列の初期化
{
	int i, j;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			Q[i][j] = 0.0;
		}
	}
	return Q[3][3];
}

double Gyoretsu_Kakesan(double Q[3][3], double R[3][3], double A[3][3])/*この行列の積（演算結果）は最後に用意された行列に入る事に注意*/
{
	int i, j, k;
	Syokika(A);

	for (i = 0; i < 3; i++)
	{
		for (k = 0; k < 3; k++)
		{
			for (j = 0; j < 3; j++)
			{
				A[i][j] += Q[i][j] * R[j][k];
			}
		}
	}
	return A[3][3];
}

double Input(double A[3][3])
{
	int i, j;
	Syokika(A);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			printf("A_%d%dを入力してください。: ", i, j);
			scanf_s("%lf", &A[i][j]);
		}
	}
	return A[3][3];
}

double Stress_Strain(double T_sigma[3][3], double Q[3][3], double T_epsilon[3][3], double result[3][3])
{
	//printf("x-y座標系の応力とひずみの関係式の係数行列\n=====================\n");
	GyakuGyoretsu(T_sigma);
	double result1[3][3];
	Syokika(result1);
	Gyoretsu_Kakesan(T_sigma, Q, result1);
	Gyoretsu_Kakesan(result1, T_epsilon, result);
	//Display(result);
	//printf("=====================\n");
	return result[3][3];
}

double In_plane_stiffness_matrix_function(double Q_bar[3][3], double Thickness, double In_plane_stiffness_matrix[3][3])
{
	/*面内剛性行列(In_plane_stiffness_matrix)を求める*/
	int i, j;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			In_plane_stiffness_matrix[i][j] += Q_bar[i][j] * Thickness;
		}
	}
	return In_plane_stiffness_matrix[3][3];
}

double Tensile_bending_coupling_matrix_function(double Q_bar[3][3], double Thickness, double Tensile_bending_coupling_matrix[3][3])
{
	int i, j;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			Tensile_bending_coupling_matrix[i][j] += Q_bar[i][j] * pow(Thickness, 2) * 0.5;
		}
	}
	return Tensile_bending_coupling_matrix[3][3];
}

double Out_of_plane_stiffness_matrix_function(double Q_bar[3][3], double Thickness, double Out_of_plane_stiffness_matrix[3][3])
{
	int i, j;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			Out_of_plane_stiffness_matrix[i][j] += Q_bar[i][j] * pow(Thickness, 3) * 0.333333;
		}
	}
	return Out_of_plane_stiffness_matrix[3][3];
}
// プログラムの実行: Ctrl + F5 または [デバッグ] > [デバッグなしで開始] メニュー
// プログラムのデバッグ: F5 または [デバッグ] > [デバッグの開始] メニュー

// 作業を開始するためのヒント: 
//    1. ソリューション エクスプローラー ウィンドウを使用してファイルを追加/管理します 
//   2. チーム エクスプローラー ウィンドウを使用してソース管理に接続します
//   3. 出力ウィンドウを使用して、ビルド出力とその他のメッセージを表示します
//   4. エラー一覧ウィンドウを使用してエラーを表示します
//   5. [プロジェクト] > [新しい項目の追加] と移動して新しいコード ファイルを作成するか、[プロジェクト] > [既存の項目の追加] と移動して既存のコード ファイルをプロジェクトに追加します
//   6. 後ほどこのプロジェクトを再び開く場合、[ファイル] > [開く] > [プロジェクト] と移動して .sln ファイルを選択します
