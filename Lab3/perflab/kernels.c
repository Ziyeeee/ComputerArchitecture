/********************************************************
 * Performance Lab: kenels.c
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"

/* 
 * Please fill in the following team struct 
 */
team_t team = {
    "王子烨",     /* 姓名 */
    "SZ170210229",  /* 学号 */
};

/***************
 * ROTATE KERNEL
 ***************/

/******************************************************
 * Your different versions of the rotate kernel go here
 ******************************************************/

/* 
 * naive_rotate - The naive baseline version of rotate 
 */
char naive_rotate_descr[] = "naive_rotate: Naive baseline implementation";
void naive_rotate(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
	    for (j = 0; j < dim; j++)
	        dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
}

/* 
 * rotate - Your current working version of rotate
 * IMPORTANT: This is the version you will be graded on
 */
void rotate_unroll_col_func_v2(int dim, pixel *src, pixel *dst);
char rotate_descr[] = "rotate: Current working version";
void rotate(int dim, pixel *src, pixel *dst) 
{
    rotate_unroll_col_func_v2(dim, src, dst);
}

/* 
 * The process of trying new methods...
 */
char rotate_descr_split[] = "rotate: 尝试分成4*4的小块，提高空间局部性";
void rotate_split(int dim, pixel *src, pixel *dst){
	int i, j, ii, jj;

    for (ii = 0; ii < dim; ii+=4)
	    for (jj = 0; jj < dim; jj+=4)
			for (i = ii; i < ii+4; i++)
				for (j = jj; j< jj+4; j++)
					dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
}

char rotate_descr_loop_unroll[] = "rotate: 在16*16分块的基础上，循环展开4*4";
void rotate_loop_unroll(int dim, pixel *src, pixel *dst){
	int i, j, ii, jj;
	// 这里16*16分块，再循环展开4*4
    for (ii = 0; ii < dim; ii+=16)
	    for (jj = 0; jj < dim; jj+=16)
			for (i = ii; i < ii+16; i+=4)
				for (j = jj; j< jj+16; j+=4){
					//相当于在4*4的小块内，手写每一个像素点的旋转变换
					dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
					dst[RIDX(dim-1-j-1, i, dim)] = src[RIDX(i, j+1, dim)];
					dst[RIDX(dim-1-j-2, i, dim)] = src[RIDX(i, j+2, dim)];
					dst[RIDX(dim-1-j-3, i, dim)] = src[RIDX(i, j+3, dim)];
					dst[RIDX(dim-1-j, i+1, dim)] = src[RIDX(i+1, j, dim)];
					dst[RIDX(dim-1-j-1, i+1, dim)] = src[RIDX(i+1, j+1, dim)];
					dst[RIDX(dim-1-j-2, i+1, dim)] = src[RIDX(i+1, j+2, dim)];
					dst[RIDX(dim-1-j-3, i+1, dim)] = src[RIDX(i+1, j+3, dim)];
					dst[RIDX(dim-1-j, i+2, dim)] = src[RIDX(i+2, j, dim)];
					dst[RIDX(dim-1-j-1, i+2, dim)] = src[RIDX(i+2, j+1, dim)];
					dst[RIDX(dim-1-j-2, i+2, dim)] = src[RIDX(i+2, j+2, dim)];
					dst[RIDX(dim-1-j-3, i+2, dim)] = src[RIDX(i+2, j+3, dim)];
					dst[RIDX(dim-1-j, i+3, dim)] = src[RIDX(i+3, j, dim)];
					dst[RIDX(dim-1-j-1, i+3, dim)] = src[RIDX(i+3, j+1, dim)];
					dst[RIDX(dim-1-j-2, i+3, dim)] = src[RIDX(i+3, j+2, dim)];
					dst[RIDX(dim-1-j-3, i+3, dim)] = src[RIDX(i+3, j+3, dim)];
				}
}


char rotate_descr_unroll_col[] = "rotate: 每16行为一块，逐列来转换";
void rotate_unroll_col(int dim, pixel *src, pixel *dst){
	int i,j;
	//每16行为一块
	for (i = 0; i < dim; i+=16)
		for (j = 0; j < dim; j++) {
			// 逐列来转换
			dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
			dst[RIDX(dim-1-j, i+1, dim)] = src[RIDX(i+1, j, dim)];
			dst[RIDX(dim-1-j, i+2, dim)] = src[RIDX(i+2, j, dim)];
			dst[RIDX(dim-1-j, i+3, dim)] = src[RIDX(i+3, j, dim)];
			dst[RIDX(dim-1-j, i+4, dim)] = src[RIDX(i+4, j, dim)];
			dst[RIDX(dim-1-j, i+5, dim)] = src[RIDX(i+5, j, dim)];
			dst[RIDX(dim-1-j, i+6, dim)] = src[RIDX(i+6, j, dim)];
			dst[RIDX(dim-1-j, i+7, dim)] = src[RIDX(i+7, j, dim)];
			dst[RIDX(dim-1-j, i+8, dim)] = src[RIDX(i+8, j, dim)];
			dst[RIDX(dim-1-j, i+9, dim)] = src[RIDX(i+9, j, dim)];
			dst[RIDX(dim-1-j, i+10, dim)] = src[RIDX(i+10, j, dim)];
			dst[RIDX(dim-1-j, i+11, dim)] = src[RIDX(i+11, j, dim)];
			dst[RIDX(dim-1-j, i+12, dim)] = src[RIDX(i+12, j, dim)];
			dst[RIDX(dim-1-j, i+13, dim)] = src[RIDX(i+13, j, dim)];
			dst[RIDX(dim-1-j, i+14, dim)] = src[RIDX(i+14, j, dim)];
			dst[RIDX(dim-1-j, i+15, dim)] = src[RIDX(i+15, j, dim)];
		}	
}


char rotate_descr_unroll_col_func[] = "rotate: 每16行为一块，逐列来转换，代入宏定义RIDX";
void rotate_unroll_col_func(int dim, pixel *src, pixel *dst){
	int i,j;
	//每16行为一块
	for (i = 0; i < dim; i+=16)
		for (j = 0; j < dim; j++) {
			// 逐列来转换
			// define RIDX(i, j, n)  ((i)*(n)+(j))
			int destmp = (dim-1-j)*dim+i;
			int srctmp = i*dim + j;
			dst[destmp] = src[srctmp]; srctmp += dim;
			dst[++destmp] = src[srctmp]; srctmp += dim;
			dst[++destmp] = src[srctmp]; srctmp += dim;
			dst[++destmp] = src[srctmp]; srctmp += dim;
			dst[++destmp] = src[srctmp]; srctmp += dim;
			dst[++destmp] = src[srctmp]; srctmp += dim;
			dst[++destmp] = src[srctmp]; srctmp += dim;
			dst[++destmp] = src[srctmp]; srctmp += dim;
			dst[++destmp] = src[srctmp]; srctmp += dim;
			dst[++destmp] = src[srctmp]; srctmp += dim;
			dst[++destmp] = src[srctmp]; srctmp += dim;
			dst[++destmp] = src[srctmp]; srctmp += dim;
			dst[++destmp] = src[srctmp]; srctmp += dim;
			dst[++destmp] = src[srctmp]; srctmp += dim;
			dst[++destmp] = src[srctmp]; srctmp += dim;
			dst[++destmp] = src[srctmp];
		}	
}

/*
char rotate_descr_unroll_row[] = "rotate: 每16列为一块，逐行来转换";
void rotate_unroll_row(int dim, pixel *src, pixel *dst){
	int i,j;
	//每16列为一块
	for (j = 0; j < dim; j+=16)
		for (i = 0; i < dim; i++) {
			//逐行来转换
			dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
			dst[RIDX(dim-2-j, i, dim)] = src[RIDX(i, j+1, dim)];
			dst[RIDX(dim-3-j, i, dim)] = src[RIDX(i, j+2, dim)];
			dst[RIDX(dim-4-j, i, dim)] = src[RIDX(i, j+3, dim)];
			dst[RIDX(dim-5-j, i, dim)] = src[RIDX(i, j+4, dim)];
			dst[RIDX(dim-6-j, i, dim)] = src[RIDX(i, j+5, dim)];
			dst[RIDX(dim-7-j, i, dim)] = src[RIDX(i, j+6, dim)];
			dst[RIDX(dim-8-j, i, dim)] = src[RIDX(i, j+7, dim)];
			dst[RIDX(dim-9-j, i, dim)] = src[RIDX(i, j+8, dim)];
			dst[RIDX(dim-10-j, i, dim)] = src[RIDX(i, j+9, dim)];
			dst[RIDX(dim-11-j, i, dim)] = src[RIDX(i, j+10, dim)];
			dst[RIDX(dim-12-j, i, dim)] = src[RIDX(i, j+11, dim)];
			dst[RIDX(dim-13-j, i, dim)] = src[RIDX(i, j+12, dim)];
			dst[RIDX(dim-14-j, i, dim)] = src[RIDX(i, j+13, dim)];
			dst[RIDX(dim-15-j, i, dim)] = src[RIDX(i, j+14, dim)];
			dst[RIDX(dim-16-j, i, dim)] = src[RIDX(i, j+15, dim)];
		}	
}
*/

char rotate_descr_unroll_col32[] = "rotate: 每32行为一块，逐列来转换";
void rotate_unroll_col32(int dim, pixel *src, pixel *dst){
	int i,j;
	//每32行为一块
	for (i = 0; i < dim; i+=32)
		for (j = 0; j < dim; j++) {
			// 逐列来转换
			dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
			dst[RIDX(dim-1-j, i+1, dim)] = src[RIDX(i+1, j, dim)];
			dst[RIDX(dim-1-j, i+2, dim)] = src[RIDX(i+2, j, dim)];
			dst[RIDX(dim-1-j, i+3, dim)] = src[RIDX(i+3, j, dim)];
			dst[RIDX(dim-1-j, i+4, dim)] = src[RIDX(i+4, j, dim)];
			dst[RIDX(dim-1-j, i+5, dim)] = src[RIDX(i+5, j, dim)];
			dst[RIDX(dim-1-j, i+6, dim)] = src[RIDX(i+6, j, dim)];
			dst[RIDX(dim-1-j, i+7, dim)] = src[RIDX(i+7, j, dim)];
			dst[RIDX(dim-1-j, i+8, dim)] = src[RIDX(i+8, j, dim)];
			dst[RIDX(dim-1-j, i+9, dim)] = src[RIDX(i+9, j, dim)];
			dst[RIDX(dim-1-j, i+10, dim)] = src[RIDX(i+10, j, dim)];
			dst[RIDX(dim-1-j, i+11, dim)] = src[RIDX(i+11, j, dim)];
			dst[RIDX(dim-1-j, i+12, dim)] = src[RIDX(i+12, j, dim)];
			dst[RIDX(dim-1-j, i+13, dim)] = src[RIDX(i+13, j, dim)];
			dst[RIDX(dim-1-j, i+14, dim)] = src[RIDX(i+14, j, dim)];
			dst[RIDX(dim-1-j, i+15, dim)] = src[RIDX(i+15, j, dim)];
			dst[RIDX(dim-1-j, i+16, dim)] = src[RIDX(i+16, j, dim)];
			dst[RIDX(dim-1-j, i+17, dim)] = src[RIDX(i+17, j, dim)];
			dst[RIDX(dim-1-j, i+18, dim)] = src[RIDX(i+18, j, dim)];
			dst[RIDX(dim-1-j, i+19, dim)] = src[RIDX(i+19, j, dim)];
			dst[RIDX(dim-1-j, i+20, dim)] = src[RIDX(i+20, j, dim)];
			dst[RIDX(dim-1-j, i+21, dim)] = src[RIDX(i+21, j, dim)];
			dst[RIDX(dim-1-j, i+22, dim)] = src[RIDX(i+22, j, dim)];
			dst[RIDX(dim-1-j, i+23, dim)] = src[RIDX(i+23, j, dim)];
			dst[RIDX(dim-1-j, i+24, dim)] = src[RIDX(i+24, j, dim)];
			dst[RIDX(dim-1-j, i+25, dim)] = src[RIDX(i+25, j, dim)];
			dst[RIDX(dim-1-j, i+26, dim)] = src[RIDX(i+26, j, dim)];
			dst[RIDX(dim-1-j, i+27, dim)] = src[RIDX(i+27, j, dim)];
			dst[RIDX(dim-1-j, i+28, dim)] = src[RIDX(i+28, j, dim)];
			dst[RIDX(dim-1-j, i+29, dim)] = src[RIDX(i+29, j, dim)];
			dst[RIDX(dim-1-j, i+30, dim)] = src[RIDX(i+30, j, dim)];
			dst[RIDX(dim-1-j, i+31, dim)] = src[RIDX(i+31, j, dim)];
		}	
}


char rotate_descr_unroll_col_func_v2[] = "rotate: 每16行为一块，逐列来转换，代入宏定义RIDX（指针实现）";
void rotate_unroll_col_func_v2(int dim, pixel *src, pixel *dst){
	int i, j;
	//src初始化为第一行第一个像素点，dst与之对应，应指向最后一行第一个像素点
	dst += (dim-1)*dim;
	for (i = 0; i < dim; i+=16){ 
		for (j = 0; j < dim; j++){ 
			//循环展开，共16次：
			*dst=*src; src+=dim; dst+=1;
			*dst=*src; src+=dim; dst+=1;
			*dst=*src; src+=dim; dst+=1;
			*dst=*src; src+=dim; dst+=1;
			*dst=*src; src+=dim; dst+=1;
			*dst=*src; src+=dim; dst+=1;
			*dst=*src; src+=dim; dst+=1;
			*dst=*src; src+=dim; dst+=1;
			*dst=*src; src+=dim; dst+=1;
			*dst=*src; src+=dim; dst+=1;
			*dst=*src; src+=dim; dst+=1;
			*dst=*src; src+=dim; dst+=1;
			*dst=*src; src+=dim; dst+=1;
			*dst=*src; src+=dim; dst+=1;
			*dst=*src; src+=dim; dst+=1;
			*dst=*src;
			//注意修改src与dst指向的位置，改为下一列的相应32行
			//int destmp = (dim-1-j)*dim+i;
			//int srctmp = i*dim + j;
			src = src-(dim<<4)+dim+1; //-15dim+1 
			dst = dst-15-dim;
		}
		dst += dim*dim+16;
		src = src+(dim<<4)-dim; //15dim
	}
}

/*********************************************************************
 * register_rotate_functions - 通过为每一个测试函数调用add_rotate_function(),
 * 登记你所有不同版本的旋转代码到评测程序driver中，
 * 当你运行driver程序时，它将测试并且给出每一个已经登记的测试函数的性能。
 *********************************************************************/

void register_rotate_functions() 
{
    add_rotate_function(&naive_rotate, naive_rotate_descr);   
    add_rotate_function(&rotate, rotate_descr);   
    /* ... Register additional test functions here */
	add_rotate_function(&rotate_split, rotate_descr_split);
	add_rotate_function(&rotate_loop_unroll, rotate_descr_loop_unroll);
	add_rotate_function(&rotate_unroll_col, rotate_descr_unroll_col);
	// add_rotate_function(&rotate_unroll_row, rotate_descr_unroll_row);//效果很差
	add_rotate_function(&rotate_unroll_col32, rotate_descr_unroll_col32);
	add_rotate_function(&rotate_unroll_col_func, rotate_descr_unroll_col_func);
	add_rotate_function(&rotate_unroll_col_func_v2, rotate_descr_unroll_col_func_v2);
}


/***************
 * SMOOTH KERNEL
 **************/

/***************************************************************
 * Various typedefs and helper functions for the smooth function
 * You may modify these any way you like.
 **************************************************************/

/* A struct used to compute averaged pixel value */
typedef struct {
    int red;
    int green;
    int blue;
    int num;
} pixel_sum;

/* Compute min and max of two integers, respectively */
static int min(int a, int b) { return (a < b ? a : b); }
static int max(int a, int b) { return (a > b ? a : b); }

/* 
 * initialize_pixel_sum - Initializes all fields of sum to 0 
 */
static void initialize_pixel_sum(pixel_sum *sum) 
{
    sum->red = sum->green = sum->blue = 0;
    sum->num = 0;
    return;
}

/* 
 * accumulate_sum - Accumulates field values of p in corresponding 
 * fields of sum 
 */
static void accumulate_sum(pixel_sum *sum, pixel p) 
{
    sum->red += (int) p.red;
    sum->green += (int) p.green;
    sum->blue += (int) p.blue;
    sum->num++;
    return;
}

/* 
 * assign_sum_to_pixel - Computes averaged pixel value in current_pixel 
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum) 
{
    current_pixel->red = (unsigned short) (sum.red/sum.num);
    current_pixel->green = (unsigned short) (sum.green/sum.num);
    current_pixel->blue = (unsigned short) (sum.blue/sum.num);
    return;
}

/* 
 * avg - Returns averaged pixel value at (i,j) 
 */
static pixel avg(int dim, int i, int j, pixel *src) 
{
    int ii, jj;
    pixel_sum sum;
    pixel current_pixel;

    initialize_pixel_sum(&sum);
    for(ii = max(i-1, 0); ii <= min(i+1, dim-1); ii++) 
	for(jj = max(j-1, 0); jj <= min(j+1, dim-1); jj++) 
	    accumulate_sum(&sum, src[RIDX(ii, jj, dim)]);

    assign_sum_to_pixel(&current_pixel, sum);
    return current_pixel;
}

/******************************************************
 * Your different versions of the smooth kernel go here
 ******************************************************/

/*
 * naive_smooth - The naive baseline version of smooth 
 */
char naive_smooth_descr[] = "naive_smooth: Naive baseline implementation";
void naive_smooth(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
	for (j = 0; j < dim; j++)
	    dst[RIDX(i, j, dim)] = avg(dim, i, j, src);
}

/*
 * smooth - Your current working version of smooth. 
 * IMPORTANT: This is the version you will be graded on
 */
void smooth_nofunc_unroll(int dim, pixel *src, pixel *dst);
char smooth_descr[] = "smooth: Current working version";
void smooth(int dim, pixel *src, pixel *dst) 
{
    smooth_nofunc_unroll(dim, src, dst);
}

char smooth_nofunc_descr[] = "smooth: 消除函数调用";
void smooth_nofunc(int dim, pixel *src, pixel *dst) 
{
    int i, j;
	//角落4个点
	//边缘6个点
	//中间8个点
	
	//左上角顶点
	dst[0].red=(src[0].red+src[1].red+src[dim].red+src[dim+1].red)/4;
	dst[0].green=(src[0].green+src[1].green+src[dim].green+src[dim+1].green)/4;
	dst[0].blue=(src[0].blue+src[1].blue+src[dim].blue+src[dim+1].blue)/4;
	
	//第一行非左上角、非右上角
	for (j = 1; j < dim-1; j++){
		dst[j].red=(src[j-1].red+src[j].red+src[j+1].red+src[dim+j-1].red+src[dim+j].red+src[dim+j+1].red)/6;
		dst[j].green=(src[j-1].green+src[j].green+src[j+1].green+src[dim+j-1].green+src[dim+j].green+src[dim+j+1].green)/6;
		dst[j].blue=(src[j-1].blue+src[j].blue+src[j+1].blue+src[dim+j-1].blue+src[dim+j].blue+src[dim+j+1].blue)/6;
	}
	
	//右上角顶点
	dst[j].red=(src[j].red+src[j-1].red+src[dim+j].red+src[dim+j-1].red)/4;
	dst[j].green=(src[j].green+src[j-1].green+src[dim+j].green+src[dim+j-1].green)/4;
	dst[j].blue=(src[j].blue+src[j-1].blue+src[dim+j].blue+src[dim+j-1].blue)/4;
	
	//逐行处理（直到倒数第2排）：左右6个，中间9个
    for (i = 1; i < dim-1; i++){
		dst[i*dim].red=(src[(i-1)*dim].red+src[(i-1)*dim+1].red+src[i*dim].red+src[i*dim+1].red+src[(i+1)*dim].red+src[(i+1)*dim+1].red)/6;
		dst[i*dim].green=(src[(i-1)*dim].green+src[(i-1)*dim+1].green+src[i*dim].green+src[i*dim+1].green+src[(i+1)*dim].green+src[(i+1)*dim+1].green)/6;
		dst[i*dim].blue=(src[(i-1)*dim].blue+src[(i-1)*dim+1].blue+src[i*dim].blue+src[i*dim+1].blue+src[(i+1)*dim].blue+src[(i+1)*dim+1].blue)/6;
		
		for (j = 1; j < dim-1; j++){
			dst[i*dim+j].red=(src[(i-1)*dim+j-1].red+src[(i-1)*dim+j].red+src[(i-1)*dim+j+1].red+src[i*dim+j-1].red+src[i*dim+j].red+src[i*dim+j+1].red+src[(i+1)*dim+j-1].red+src[(i+1)*dim+j].red+src[(i+1)*dim+j+1].red)/9;
			dst[i*dim+j].green=(src[(i-1)*dim+j-1].green+src[(i-1)*dim+j].green+src[(i-1)*dim+j+1].green+src[i*dim+j-1].green+src[i*dim+j].green+src[i*dim+j+1].green+src[(i+1)*dim+j-1].green+src[(i+1)*dim+j].green+src[(i+1)*dim+j+1].green)/9;
			dst[i*dim+j].blue=(src[(i-1)*dim+j-1].blue+src[(i-1)*dim+j].blue+src[(i-1)*dim+j+1].blue+src[i*dim+j-1].blue+src[i*dim+j].blue+src[i*dim+j+1].blue+src[(i+1)*dim+j-1].blue+src[(i+1)*dim+j].blue+src[(i+1)*dim+j+1].blue)/9;
		}
		
		dst[i*dim+j].red=(src[(i-1)*dim+j-1].red+src[(i-1)*dim+j].red+src[i*dim+j-1].red+src[i*dim+j].red+src[(i+1)*dim+j-1].red+src[(i+1)*dim+j].red)/6;
		dst[i*dim+j].green=(src[(i-1)*dim+j-1].green+src[(i-1)*dim+j].green+src[i*dim+j-1].green+src[i*dim+j].green+src[(i+1)*dim+j-1].green+src[(i+1)*dim+j].green)/6;
		dst[i*dim+j].blue=(src[(i-1)*dim+j-1].blue+src[(i-1)*dim+j].blue+src[i*dim+j-1].blue+src[i*dim+j].blue+src[(i+1)*dim+j-1].blue+src[(i+1)*dim+j].blue)/6;
	}
	
	//左下角顶点
	dst[i*dim].red=(src[(i-1)*dim].red+src[(i-1)*dim+1].red+src[i*dim].red+src[i*dim+1].red)/4;
	dst[i*dim].green=(src[(i-1)*dim].green+src[(i-1)*dim+1].green+src[i*dim].green+src[i*dim+1].green)/4;
	dst[i*dim].blue=(src[(i-1)*dim].blue+src[(i-1)*dim+1].blue+src[i*dim].blue+src[i*dim+1].blue)/4;
	
	//最后一行非左下角、非右下角
	for (j = 1; j < dim-1; j++) {
		dst[i*dim+j].red=(src[(i-1)*dim+j-1].red+src[(i-1)*dim+j].red+src[(i-1)*dim+j+1].red+src[i*dim+j-1].red+src[i*dim+j].red+src[i*dim+j+1].red)/6;
		dst[i*dim+j].green=(src[(i-1)*dim+j-1].green+src[(i-1)*dim+j].green+src[(i-1)*dim+j+1].green+src[i*dim+j-1].green+src[i*dim+j].green+src[i*dim+j+1].green)/6;
		dst[i*dim+j].blue=(src[(i-1)*dim+j-1].blue+src[(i-1)*dim+j].blue+src[(i-1)*dim+j+1].blue+src[i*dim+j-1].blue+src[i*dim+j].blue+src[i*dim+j+1].blue)/6;
	}
	
	//右下角顶点
	dst[i*dim+j].red=(src[(i-1)*dim+j-1].red+src[(i-1)*dim+j].red+src[i*dim+j-1].red+src[i*dim+j].red)/4;
	dst[i*dim+j].green=(src[(i-1)*dim+j-1].green+src[(i-1)*dim+j].green+src[i*dim+j-1].green+src[i*dim+j].green)/4;
	dst[i*dim+j].blue=(src[(i-1)*dim+j-1].blue+src[(i-1)*dim+j].blue+src[i*dim+j-1].blue+src[i*dim+j].blue)/4;	
}


char smooth_nofunc_pointer_descr[] = "smooth: 消除函数调用+减少重复计算";
void smooth_nofunc_pointer(int dim, pixel *src, pixel *dst){
	int i, j;
	//以下处理4个角落+4条边
	dst[0].red=(src[0].red+src[1].red+src[dim].red+src[dim+1].red)/4;
	dst[0].green=(src[0].green+src[1].green+src[dim].green+src[dim+1].green)/4;
	dst[0].blue=(src[0].blue+src[1].blue+src[dim].blue+src[dim+1].blue)/4;
	
	for (j = 1; j < dim-1; j++){
		dst[j].red=(src[j-1].red+src[j].red+src[j+1].red+src[dim+j-1].red+src[dim+j].red+src[dim+j+1].red)/6;
		dst[j].green=(src[j-1].green+src[j].green+src[j+1].green+src[dim+j-1].green+src[dim+j].green+src[dim+j+1].green)/6;
		dst[j].blue=(src[j-1].blue+src[j].blue+src[j+1].blue+src[dim+j-1].blue+src[dim+j].blue+src[dim+j+1].blue)/6;
	}
	
	dst[j].red=(src[j].red+src[j-1].red+src[dim+j].red+src[dim+j-1].red)/4;
	dst[j].green=(src[j].green+src[j-1].green+src[dim+j].green+src[dim+j-1].green)/4;
	dst[j].blue=(src[j].blue+src[j-1].blue+src[dim+j].blue+src[dim+j-1].blue)/4;
	
    for (i = 1; i < dim-1; i++){
		dst[i*dim].red=(src[(i-1)*dim].red+src[(i-1)*dim+1].red+src[i*dim].red+src[i*dim+1].red+src[(i+1)*dim].red+src[(i+1)*dim+1].red)/6;
		dst[i*dim].green=(src[(i-1)*dim].green+src[(i-1)*dim+1].green+src[i*dim].green+src[i*dim+1].green+src[(i+1)*dim].green+src[(i+1)*dim+1].green)/6;
		dst[i*dim].blue=(src[(i-1)*dim].blue+src[(i-1)*dim+1].blue+src[i*dim].blue+src[i*dim+1].blue+src[(i+1)*dim].blue+src[(i+1)*dim+1].blue)/6;
		
		dst[i*dim+j].red=(src[(i-1)*dim+j-1].red+src[(i-1)*dim+j].red+src[i*dim+j-1].red+src[i*dim+j].red+src[(i+1)*dim+j-1].red+src[(i+1)*dim+j].red)/6;
		dst[i*dim+j].green=(src[(i-1)*dim+j-1].green+src[(i-1)*dim+j].green+src[i*dim+j-1].green+src[i*dim+j].green+src[(i+1)*dim+j-1].green+src[(i+1)*dim+j].green)/6;
		dst[i*dim+j].blue=(src[(i-1)*dim+j-1].blue+src[(i-1)*dim+j].blue+src[i*dim+j-1].blue+src[i*dim+j].blue+src[(i+1)*dim+j-1].blue+src[(i+1)*dim+j].blue)/6;
	}
	
	dst[i*dim].red=(src[(i-1)*dim].red+src[(i-1)*dim+1].red+src[i*dim].red+src[i*dim+1].red)/4;
	dst[i*dim].green=(src[(i-1)*dim].green+src[(i-1)*dim+1].green+src[i*dim].green+src[i*dim+1].green)/4;
	dst[i*dim].blue=(src[(i-1)*dim].blue+src[(i-1)*dim+1].blue+src[i*dim].blue+src[i*dim+1].blue)/4;
	
	for (j = 1; j < dim-1; j++) {
		dst[i*dim+j].red=(src[(i-1)*dim+j-1].red+src[(i-1)*dim+j].red+src[(i-1)*dim+j+1].red+src[i*dim+j-1].red+src[i*dim+j].red+src[i*dim+j+1].red)/6;
		dst[i*dim+j].green=(src[(i-1)*dim+j-1].green+src[(i-1)*dim+j].green+src[(i-1)*dim+j+1].green+src[i*dim+j-1].green+src[i*dim+j].green+src[i*dim+j+1].green)/6;
		dst[i*dim+j].blue=(src[(i-1)*dim+j-1].blue+src[(i-1)*dim+j].blue+src[(i-1)*dim+j+1].blue+src[i*dim+j-1].blue+src[i*dim+j].blue+src[i*dim+j+1].blue)/6;
	}
	
	dst[i*dim+j].red=(src[(i-1)*dim+j-1].red+src[(i-1)*dim+j].red+src[i*dim+j-1].red+src[i*dim+j].red)/4;
	dst[i*dim+j].green=(src[(i-1)*dim+j-1].green+src[(i-1)*dim+j].green+src[i*dim+j-1].green+src[i*dim+j].green)/4;
	dst[i*dim+j].blue=(src[(i-1)*dim+j-1].blue+src[(i-1)*dim+j].blue+src[i*dim+j-1].blue+src[i*dim+j].blue)/4;	
	
	//下面减少重复计算+用指针实现
	pixel *row1 = &src[0];
	pixel *row2 = &src[dim];
	pixel *row3 = &src[dim+dim];
	dst = &dst[dim+1];
	int left_red, mid_red, right_red;
	int left_green, mid_green, right_green;
	int left_blue, mid_blue, right_blue;
	for (i=1; i<dim-1; i++){
		left_red = row1->red + row2->red + row3->red;
		left_green = row1->green + row2->green + row3->green;
		left_blue = row1->blue + row2->blue + row3->blue;
		row1++;
		row2++;
		row3++;
		
		mid_red = row1->red + row2->red + row3->red;
		mid_green = row1->green + row2->green + row3->green;
		mid_blue = row1->blue + row2->blue + row3->blue;
		row1++;
		row2++;
		row3++;
		
		for (j=1; j<dim-1; j++){
			right_red = row1->red + row2->red + row3->red;
			right_green = row1->green + row2->green + row3->green;
			right_blue = row1->blue + row2->blue + row3->blue;
			row1++;
			row2++;
			row3++;
			
			dst->red = (left_red + mid_red + right_red)/9;
			dst->green = (left_green + mid_green + right_green)/9;
			dst->blue = (left_blue + mid_blue + right_blue)/9;
			dst++;
			
			left_red = mid_red;
			left_blue = mid_blue;
			left_green = mid_green;
			
			mid_red = right_red;
			mid_blue = right_blue;
			mid_green = right_green;
			
		}
		dst += 2;
	}
	
}

char smooth_nofunc_unroll_descr[] = "smooth: 消除函数调用+减少重复计算+循环展开";
void smooth_nofunc_unroll(int dim, pixel *src, pixel *dst){
	int i, j;
	//以下处理4个角落+4条边
	dst[0].red=(src[0].red+src[1].red+src[dim].red+src[dim+1].red)/4;
	dst[0].green=(src[0].green+src[1].green+src[dim].green+src[dim+1].green)/4;
	dst[0].blue=(src[0].blue+src[1].blue+src[dim].blue+src[dim+1].blue)/4;
	
	for (j = 1; j < dim-1; j++){
		dst[j].red=(src[j-1].red+src[j].red+src[j+1].red+src[dim+j-1].red+src[dim+j].red+src[dim+j+1].red)/6;
		dst[j].green=(src[j-1].green+src[j].green+src[j+1].green+src[dim+j-1].green+src[dim+j].green+src[dim+j+1].green)/6;
		dst[j].blue=(src[j-1].blue+src[j].blue+src[j+1].blue+src[dim+j-1].blue+src[dim+j].blue+src[dim+j+1].blue)/6;
	}
	
	dst[j].red=(src[j].red+src[j-1].red+src[dim+j].red+src[dim+j-1].red)/4;
	dst[j].green=(src[j].green+src[j-1].green+src[dim+j].green+src[dim+j-1].green)/4;
	dst[j].blue=(src[j].blue+src[j-1].blue+src[dim+j].blue+src[dim+j-1].blue)/4;
	
    for (i = 1; i < dim-1; i++){
		dst[i*dim].red=(src[(i-1)*dim].red+src[(i-1)*dim+1].red+src[i*dim].red+src[i*dim+1].red+src[(i+1)*dim].red+src[(i+1)*dim+1].red)/6;
		dst[i*dim].green=(src[(i-1)*dim].green+src[(i-1)*dim+1].green+src[i*dim].green+src[i*dim+1].green+src[(i+1)*dim].green+src[(i+1)*dim+1].green)/6;
		dst[i*dim].blue=(src[(i-1)*dim].blue+src[(i-1)*dim+1].blue+src[i*dim].blue+src[i*dim+1].blue+src[(i+1)*dim].blue+src[(i+1)*dim+1].blue)/6;
		
		dst[i*dim+j].red=(src[(i-1)*dim+j-1].red+src[(i-1)*dim+j].red+src[i*dim+j-1].red+src[i*dim+j].red+src[(i+1)*dim+j-1].red+src[(i+1)*dim+j].red)/6;
		dst[i*dim+j].green=(src[(i-1)*dim+j-1].green+src[(i-1)*dim+j].green+src[i*dim+j-1].green+src[i*dim+j].green+src[(i+1)*dim+j-1].green+src[(i+1)*dim+j].green)/6;
		dst[i*dim+j].blue=(src[(i-1)*dim+j-1].blue+src[(i-1)*dim+j].blue+src[i*dim+j-1].blue+src[i*dim+j].blue+src[(i+1)*dim+j-1].blue+src[(i+1)*dim+j].blue)/6;
	}
	
	dst[i*dim].red=(src[(i-1)*dim].red+src[(i-1)*dim+1].red+src[i*dim].red+src[i*dim+1].red)/4;
	dst[i*dim].green=(src[(i-1)*dim].green+src[(i-1)*dim+1].green+src[i*dim].green+src[i*dim+1].green)/4;
	dst[i*dim].blue=(src[(i-1)*dim].blue+src[(i-1)*dim+1].blue+src[i*dim].blue+src[i*dim+1].blue)/4;
	
	for (j = 1; j < dim-1; j++) {
		dst[i*dim+j].red=(src[(i-1)*dim+j-1].red+src[(i-1)*dim+j].red+src[(i-1)*dim+j+1].red+src[i*dim+j-1].red+src[i*dim+j].red+src[i*dim+j+1].red)/6;
		dst[i*dim+j].green=(src[(i-1)*dim+j-1].green+src[(i-1)*dim+j].green+src[(i-1)*dim+j+1].green+src[i*dim+j-1].green+src[i*dim+j].green+src[i*dim+j+1].green)/6;
		dst[i*dim+j].blue=(src[(i-1)*dim+j-1].blue+src[(i-1)*dim+j].blue+src[(i-1)*dim+j+1].blue+src[i*dim+j-1].blue+src[i*dim+j].blue+src[i*dim+j+1].blue)/6;
	}
	
	dst[i*dim+j].red=(src[(i-1)*dim+j-1].red+src[(i-1)*dim+j].red+src[i*dim+j-1].red+src[i*dim+j].red)/4;
	dst[i*dim+j].green=(src[(i-1)*dim+j-1].green+src[(i-1)*dim+j].green+src[i*dim+j-1].green+src[i*dim+j].green)/4;
	dst[i*dim+j].blue=(src[(i-1)*dim+j-1].blue+src[(i-1)*dim+j].blue+src[i*dim+j-1].blue+src[i*dim+j].blue)/4;	
	
	//下面用减少重复计算+指针实现
	pixel *row1 = &src[0];
	pixel *row2 = &src[dim];
	pixel *row3 = &src[dim+dim];
	dst = &dst[dim+1];
	int left_red, mid_red, right_red;
	int left_green, mid_green, right_green;
	int left_blue, mid_blue, right_blue;
	
	for (i=1; i<dim-1; i++){
		left_red = row1->red + row2->red + row3->red;
		left_green = row1->green + row2->green + row3->green;
		left_blue = row1->blue + row2->blue + row3->blue;
		row1++;
		row2++;
		row3++;
		
		mid_red = row1->red + row2->red + row3->red;
		mid_green = row1->green + row2->green + row3->green;
		mid_blue = row1->blue + row2->blue + row3->blue;
		row1++;
		row2++;
		row3++;
		// 这里进行循环展开
		for (j=1; j<dim-1; j+=2){
			right_red = row1->red + row2->red + row3->red;
			right_green = row1->green + row2->green + row3->green;
			right_blue = row1->blue + row2->blue + row3->blue;
			row1++;
			row2++;
			row3++;
			
			dst->red = (left_red + mid_red + right_red)/9;
			dst->green = (left_green + mid_green + right_green)/9;
			dst->blue = (left_blue + mid_blue + right_blue)/9;
			dst++;
			
			left_red = mid_red;
			left_blue = mid_blue;
			left_green = mid_green;
			
			mid_red = right_red;
			mid_blue = right_blue;
			mid_green = right_green;
			
			right_red = row1->red + row2->red + row3->red;
			right_green = row1->green + row2->green + row3->green;
			right_blue = row1->blue + row2->blue + row3->blue;
			row1++;
			row2++;
			row3++;
			
			dst->red = (left_red + mid_red + right_red)/9;
			dst->green = (left_green + mid_green + right_green)/9;
			dst->blue = (left_blue + mid_blue + right_blue)/9;
			dst++;
			
			left_red = mid_red;
			left_blue = mid_blue;
			left_green = mid_green;
			
			mid_red = right_red;
			mid_blue = right_blue;
			mid_green = right_green;
			
		}
		dst += 2;
	}
	
}




/********************************************************************* 
 * register_smooth_functions - 通过为每一个测试函数调用add_smooth_funtion(),
 * 登记所有不同版本的smooth代码到评测程序driver中。
 * 当你运行driver程序时，它将测试并且给出每一个已经登记的测试函数的性能。
 *********************************************************************/

void register_smooth_functions() {
	add_smooth_function(&naive_smooth, naive_smooth_descr);
    add_smooth_function(&smooth, smooth_descr);
    /* ... Register additional test functions here */
	add_smooth_function(&smooth_nofunc, smooth_nofunc_descr);
	add_smooth_function(&smooth_nofunc_pointer, smooth_nofunc_pointer_descr);
	add_smooth_function(&smooth_nofunc_unroll, smooth_nofunc_unroll_descr);
}

