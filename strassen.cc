#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int cross_point = 128;

int** create_mat(int dim) {
	int** mat = 0;
	mat = (int**) malloc(sizeof(int*) * dim);
	for(int i = 0; i < dim; i++) {
		mat[i] = (int*) malloc(sizeof(int) * dim);
		for(int j = 0; j < dim; j++) {
			mat[i][j] = 0;
		}
	}
	return mat;
}

void free_mat(int dim, int** mat) {
	for(int i = 0; i < dim; i++) {
		free(mat[i]);
	}
	free(mat);
}

int** matrix_mult(int dim, int** matrix_1, int** matrix_2) {
	int** prod_mat = create_mat(dim);
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j< dim; j++) {
			for (int k = 0; k < dim; k++) {
				prod_mat[i][j] += matrix_1[i][k]*matrix_2[k][j];
			}
		}
	}
	return prod_mat;
}

int** matrix_add(int dim, int add_sub, int** matrix_1, int** matrix_2) {
	int** sum_mat = create_mat(dim);
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			sum_mat[i][j] = matrix_1[i][j] + add_sub*matrix_2[i][j];
		}
	}
	return sum_mat;
}

int** matrix_combine(int dim, int** top_l, int** top_r, int** bot_l, int** bot_r) {
	int** mat = create_mat(2*dim);
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			mat[i][j] = top_l[i][j];
		}
	}
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			mat[i+dim][j] = bot_l[i][j];
		}
	}
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			mat[i][j+dim] = top_r[i][j];
		}
	}
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			mat[i+dim][j+dim] = bot_r[i][j];
		}
	}
	return mat;
}

int** mat_split(int dim, int quadrant, int** mat) {
	int new_dim = dim;
	if (new_dim % 2 != 0) {
		new_dim += 1;
	}
	int** split_mat = create_mat(new_dim);
	if (quadrant == 1) {
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				split_mat[i][j] = mat[i][j];
			}
		}
	}
	else if (quadrant == 2) {
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				split_mat[i][j] = mat[i][j+dim];
			}
		}

	}
	else if (quadrant == 3) {
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				split_mat[i][j] = mat[i+dim][j];
			}
		}

	}
	else {
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				split_mat[i][j] = mat[i+dim][j+dim];
			}
		}
	}
	return split_mat;
}

int** strassen(int dim, int** matrix_1, int** matrix_2) {
	if (dim <= cross_point) {
		return matrix_mult(dim, matrix_1, matrix_2);
	}
	else {
		int new_dim = dim/2;
		int fix_new_dim = new_dim;
		if (fix_new_dim % 2 != 0) {
			fix_new_dim += 1;
		}

		int** A = mat_split(new_dim, 1, matrix_1);
		int** B = mat_split(new_dim, 2, matrix_1);
		int** C = mat_split(new_dim, 3, matrix_1);
		int** D = mat_split(new_dim, 4, matrix_1);

		int** E = mat_split(new_dim, 1, matrix_2);
		int** F = mat_split(new_dim, 2, matrix_2);
		int** G = mat_split(new_dim, 3, matrix_2);
		int** H = mat_split(new_dim, 4, matrix_2);

		int** P1 = strassen(fix_new_dim, A, matrix_add(fix_new_dim, -1, F, H));
		int** P2 = strassen(fix_new_dim, matrix_add(fix_new_dim, 1, A, B), H);
		int** P3 = strassen(fix_new_dim, matrix_add(fix_new_dim, 1, C, D), E);
		int** P4 = strassen(fix_new_dim, D, matrix_add(fix_new_dim, -1, G, E));
		int** P5 = strassen(fix_new_dim, matrix_add(fix_new_dim, 1, A, D), matrix_add(fix_new_dim, 1, E, H));
		int** P6 = strassen(fix_new_dim, matrix_add(fix_new_dim, -1, B, D), matrix_add(fix_new_dim, 1, G, H));
		int** P7 = strassen(fix_new_dim, matrix_add(fix_new_dim, -1, A, C), matrix_add(fix_new_dim, 1, E, F));

		int** top_l = matrix_add(fix_new_dim, 1, P5, matrix_add(fix_new_dim, 1, P6, matrix_add(fix_new_dim, -1, P4, P2))); 
		int** top_r = matrix_add(fix_new_dim, 1, P1, P2);
		int** bot_l = matrix_add(fix_new_dim, 1, P3, P4);
		int** bot_r = matrix_add(fix_new_dim, 1, P5, matrix_add(fix_new_dim, -1, P1, matrix_add(fix_new_dim, 1, P3, P7)));

		return matrix_combine(new_dim, top_l, top_r, bot_l, bot_r);
	}
} 


int main(int argc, char *argv[]) {
	if (argc != 4) {
		fprintf(stderr, "Error: expecting 3 arguments");
		return 1;
	}
	int dim = atoi(argv[2]);
	int new_dim = dim;
	FILE *file = fopen(argv[3], "r");
	if (file == 0) {
		fprintf(stderr, "Could not open input file");
		return 1;
	}
	clock_t begin = clock();
	if (new_dim % 2 != 0) {
		new_dim += 1;
	}

	int** matrix_1 = create_mat(new_dim);
	int** matrix_2 = create_mat(new_dim);
	
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			fscanf(file, "%i", &matrix_1[i][j]);
		}
	}
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			fscanf(file, "%i", &matrix_2[i][j]);
		}
	}
	int** result_mat = strassen(dim, matrix_1, matrix_2);
	for (int i = 0; i < dim; i++) {
		fprintf(stdout, "%d\n", result_mat[i][i]);
	}
	

	clock_t end = clock();
	double time = (double)(end - begin) / CLOCKS_PER_SEC;
	// fprintf(stdout, "time taken %f", time);
}


