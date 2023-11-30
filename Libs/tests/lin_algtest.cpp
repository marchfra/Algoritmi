#include "test_config.hpp"
#include "../include/lin_alg.hpp"

TEST_CASE("testing swapRowsLinSystem function") {
	static const int N = 4;
	double **M;
	M = new double*[N];
	M[0] = new double[N*N];
	for (int i = 1; i < N; i++) {
		M[i] = M[i - 1] + N;
	}

	SUBCASE("3x3 matrix") {
		M[0][0] = 2;	M[0][1] = -1;	M[0][2] = 0;
		M[1][0] = -1;	M[1][1] = 2;	M[1][2] = -1;
		M[2][0] = 0;	M[2][1] = -1;	M[2][2] = 2;
		double v[] = {1, 2, 1};

		swapRowsLinSystem(M, v, N, 0, 1);

		CHECK(M[0][0] == -1);	CHECK(M[0][1] == 2);	CHECK(M[0][2] == -1);
		CHECK(M[1][0] == 2);	CHECK(M[1][1] == -1);	CHECK(M[1][2] == 0);
		CHECK(M[2][0] == 0);	CHECK(M[2][1] == -1);	CHECK(M[2][2] == 2);
	}

	SUBCASE("4x4 matrix") {
		M[0][0] = 1;	M[0][1] = 2;	M[0][2] = 1;	M[0][3] = -1;
		M[1][0] = 3;	M[1][1] = 2;	M[1][2] = 4;	M[1][3] = 4;
		M[2][0] = 4;	M[2][1] = 4;	M[2][2] = 3;	M[2][3] = 4;
		M[3][0] = 2;	M[3][1] = 0;	M[3][2] = 1;	M[3][3] = 5;
		double v[] = {5, 16, 22, 15};

		swapRowsLinSystem(M, v, N, 1, 3);

		CHECK(M[0][0] == 1);	CHECK(M[0][1] == 2);	CHECK(M[0][2] == 1);	CHECK(M[0][3] == -1);
		CHECK(M[1][0] == 2);	CHECK(M[1][1] == 0);	CHECK(M[1][2] == 1);	CHECK(M[1][3] == 5);
		CHECK(M[2][0] == 4);	CHECK(M[2][1] == 4);	CHECK(M[2][2] == 3);	CHECK(M[2][3] == 4);
		CHECK(M[3][0] == 3);	CHECK(M[3][1] == 2);	CHECK(M[3][2] == 4);	CHECK(M[3][3] == 4);
	}

	delete[] M[0];
	delete[] M;
}

TEST_CASE("testing partialPivoting function") {
	static const int N = 4;
	double **M;
	M = new double*[N];
	M[0] = new double[N*N];
	for (int i = 1; i < N; i++) {
		M[i] = M[i - 1] + N;
	}

	SUBCASE("3x3 matrix") {
		M[0][0] = 2;	M[0][1] = -1;	M[0][2] = 0;
		M[1][0] = -1;	M[1][1] = 2;	M[1][2] = -1;
		M[2][0] = 0;	M[2][1] = -4;	M[2][2] = 2;
		double v[] = {1, 2, 1};

		partialPivoting(M, v, 3, 1);

		CHECK(M[0][0] == 2);	CHECK(M[0][1] == -1);	CHECK(M[0][2] == 0);
		CHECK(M[1][0] == 0);	CHECK(M[1][1] == -4);	CHECK(M[1][2] == 2);
		CHECK(M[2][0] == -1);	CHECK(M[2][1] == 2);	CHECK(M[2][2] == -1);
	}

	SUBCASE("4x4 matrix") {
		M[0][0] = 1;	M[0][1] = 2;	M[0][2] = 1;	M[0][3] = -1;
		M[1][0] = 3;	M[1][1] = 2;	M[1][2] = 4;	M[1][3] = 4;
		M[2][0] = 4;	M[2][1] = 4;	M[2][2] = 3;	M[2][3] = 4;
		M[3][0] = 2;	M[3][1] = 0;	M[3][2] = 1;	M[3][3] = 5;
		double v[] = {5, 16, 22, 15};

		partialPivoting(M, v, N, 0);

		CHECK(M[0][0] == 4);	CHECK(M[0][1] == 4);	CHECK(M[0][2] == 3);	CHECK(M[0][3] == 4);
		CHECK(M[1][0] == 3);	CHECK(M[1][1] == 2);	CHECK(M[1][2] == 4);	CHECK(M[1][3] == 4);
		CHECK(M[2][0] == 1);	CHECK(M[2][1] == 2);	CHECK(M[2][2] == 1);	CHECK(M[2][3] == -1);
		CHECK(M[3][0] == 2);	CHECK(M[3][1] == 0);	CHECK(M[3][2] == 1);	CHECK(M[3][3] == 5);
	}

	delete[] M[0];
	delete[] M;
}

TEST_CASE("testing gaussElim function") {
	static const int N = 4;
	double **M;
	M = new double*[N];
	M[0] = new double[N*N];
	for (int i = 1; i < N; i++) {
		M[i] = M[i - 1] + N;
	}

	SUBCASE("3x3 matrix") {
		M[0][0] = 2;	M[0][1] = -1;	M[0][2] = 0;
		M[1][0] = -1;	M[1][1] = 2;	M[1][2] = -1;
		M[2][0] = 0;	M[2][1] = -1;	M[2][2] = 2;
		double v[] = {1, 2, 1};

		gaussElim(M, v, 3);

		CHECK(M[0][0] == doctest::Approx(2));
		CHECK(M[0][1] == doctest::Approx(-1));
		CHECK(M[0][2] == doctest::Approx(0));

		CHECK(M[1][0] == doctest::Approx(0));
		CHECK(M[1][1] == doctest::Approx(1.5));
		CHECK(M[1][2] == doctest::Approx(-1));

		CHECK(M[2][0] == doctest::Approx(0));
		CHECK(M[2][1] == doctest::Approx(0));
		CHECK(M[2][2] == doctest::Approx(4.0 / 3.0));
	}

	SUBCASE("4x4 matrix") {
		M[0][0] = 1;	M[0][1] = 2;	M[0][2] = 1;	M[0][3] = -1;
		M[1][0] = 3;	M[1][1] = 2;	M[1][2] = 4;	M[1][3] = 4;
		M[2][0] = 4;	M[2][1] = 4;	M[2][2] = 3;	M[2][3] = 4;
		M[3][0] = 2;	M[3][1] = 0;	M[3][2] = 1;	M[3][3] = 5;
		double v[] = {5, 16, 22, 15};

		gaussElim(M, v, N);

		CHECK(M[0][0] == doctest::Approx(4));
		CHECK(M[0][1] == doctest::Approx(4));
		CHECK(M[0][2] == doctest::Approx(3));
		CHECK(M[0][3] == doctest::Approx(4));

		CHECK(M[1][0] == doctest::Approx(0));
		CHECK(M[1][1] == doctest::Approx(-2));
		CHECK(M[1][2] == doctest::Approx(-0.5));
		CHECK(M[1][3] == doctest::Approx(3));

		CHECK(M[2][0] == doctest::Approx(0));
		CHECK(M[2][1] == doctest::Approx(0));
		CHECK(M[2][2] == doctest::Approx(2));
		CHECK(M[2][3] == doctest::Approx(-0.5));

		CHECK(M[3][0] == doctest::Approx(0));
		CHECK(M[3][1] == doctest::Approx(0));
		CHECK(M[3][2] == doctest::Approx(0));
		CHECK(M[3][3] == doctest::Approx(-0.5));
	}

	delete[] M[0];
	delete[] M;
}

TEST_CASE("testing solveLinSystem function") {
	static const int N = 4;
	double **M;
	M = new double*[N];
	M[0] = new double[N*N];
	for (int i = 1; i < N; i++) {
		M[i] = M[i - 1] + N;
	}
	double x[N];

	SUBCASE("3x3 matrix") {
		M[0][0] = 2;	M[0][1] = -1;	M[0][2] = 0;
		M[1][0] = -1;	M[1][1] = 2;	M[1][2] = -1;
		M[2][0] = 0;	M[2][1] = -1;	M[2][2] = 2;
		double v[] = {1, 2, 1};

		double solution[] = {2, 3, 2};

		solveLinSystem(M, v, x, 3);

		for (int i = 0; i < 3; i++) {
			CHECK(x[i] == doctest::Approx(solution[i]));
		}
	}

	SUBCASE("4x4 matrix") {
		M[0][0] = 1;	M[0][1] = 2;	M[0][2] = 1;	M[0][3] = -1;
		M[1][0] = 3;	M[1][1] = 2;	M[1][2] = 4;	M[1][3] = 4;
		M[2][0] = 4;	M[2][1] = 4;	M[2][2] = 3;	M[2][3] = 4;
		M[3][0] = 2;	M[3][1] = 0;	M[3][2] = 1;	M[3][3] = 5;
		double v[] = {5, 16, 22, 15};

		double solution[] = {16, -6, -2, -3};

		solveLinSystem(M, v, x, N);

		for (int i = 0; i < N; i++) {
			CHECK(x[i] == doctest::Approx(solution[i]));
		}
	}

	SUBCASE("4x4 matrix with /0") {
		M[0][0] = 1;	M[0][1] = 2;	M[0][2] = 1;	M[0][3] = -1;
		M[1][0] = 3;	M[1][1] = 6;	M[1][2] = 4;	M[1][3] = 4;
		M[2][0] = 4;	M[2][1] = 4;	M[2][2] = 3;	M[2][3] = 4;
		M[3][0] = 2;	M[3][1] = 0;	M[3][2] = 1;	M[3][3] = 5;
		double v[] = {5, 16, 22, 15};

		double solution[] = {4, -12, 22, -3};

		solveLinSystem(M, v, x, N);

		for (int i = 0; i < N; i++) {
			CHECK(x[i] == doctest::Approx(solution[i]));
		}
	}

	delete[] M[0];
	delete[] M;
}

TEST_CASE("testing tridiagonalSolver function") {
	const int nEq = 5;
	double a[nEq] = {nan(""), 1, 1, 1, 1};
	double b[nEq] = {2, 2, 2, 2, 2};
	double c[nEq] = {1, 1, 1, 1, nan("")};
	double v[nEq] = {1, 0, 3, 1, 0};
	double x[nEq];

	SUBCASE ("everything working properly") {
		double expected[nEq] = {2, -3, 4, -2, 1};
		tridiagonalSolver(a, b, c, v, x, nEq);

		for (int i = 0; i < nEq; i++) {
			CHECK(x[i] == doctest::Approx(expected[i]));
		}
	}

	SUBCASE("exceptions") {
		CHECK_THROWS_WITH_AS(tridiagonalSolver(a, b, c, v, x, 5000),
							 "Number of equations must not be greater than 4096.",
							 std::invalid_argument);

		double A[] = {0, 1, 1, 1, 1};
		CHECK_THROWS_WITH_AS(tridiagonalSolver(A, b, c, v, x, nEq),
							 "The first element of a must be nan(\"\").",
							 std::invalid_argument);

		double C[] = {1, 1, 1, 1, 0};
		CHECK_THROWS_WITH_AS(tridiagonalSolver(a, b, C, v, x, nEq),
							 "The last element of c must be nan(\"\").",
							 std::invalid_argument);
	}
}
