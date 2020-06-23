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
char rotate_descr[] = "rotate: Current working version";
void rotate(int dim, pixel *src, pixel *dst) 
{
    int i, j, ii, jj1, jj2;
    int block=16;//blocking the Matrix
    for(i=0; i<dim; i+=block)
    {
        for(j=0; j<dim; j+=block)
        {
            //block*block mini matrix
            for(ii=i; ii<i+block; ii++) 
            {
                for(jj1=j; jj1<j+block; jj1+=2)
                {
                    dst[RIDX(dim-1-jj1, ii, dim)] = src[RIDX(ii, jj1, dim)];
                }
                for(jj2=j+1; jj2<j+block; jj2+=2)
                {
                    dst[RIDX(dim-1-jj2, ii, dim)] = src[RIDX(ii, jj2, dim)];
                }
            }
        }
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
char smooth_descr[] = "smooth: Current working version";
void smooth(int dim, pixel *src, pixel *dst) 
{
    int i, j, pos;
    int ii, jj;
    pixel sub_sum[3];

    // left_top
    dst[0].red = (src[0].red + src[1].red + src[dim].red + src[dim+1].red) >> 2;
    dst[0].green = (src[0].green + src[1].green + src[dim].green + src[dim+1].green) >> 2;
    dst[0].blue = (src[0].blue + src[1].blue + src[dim].blue + src[dim+1].blue) >> 2;

    // right_top
    dst[dim-1].red = (src[dim-2].red + src[dim-1].red + src[dim*2-2].red + src[dim*2-1].red) >> 2;
    dst[dim-1].green = (src[dim-2].green + src[dim-1].green + src[dim*2-2].green + src[dim*2-1].green) >> 2;
    dst[dim-1].blue = (src[dim-2].blue + src[dim-1].blue + src[dim*2-2].blue + src[dim*2-1].blue) >> 2;

    // left_bottom
    dst[dim*dim-dim].red = (src[dim*dim-dim-dim].red + src[dim*dim-dim-dim+1].red + src[dim*dim-dim].red + src[dim*dim-dim+1].red) >> 2;
    dst[dim*dim-dim].green = (src[dim*dim-dim-dim].green + src[dim*dim-dim-dim+1].green + src[dim*dim-dim].green + src[dim*dim-dim+1].green) >> 2;
    dst[dim*dim-dim].blue = (src[dim*dim-dim-dim].blue + src[dim*dim-dim-dim+1].blue + src[dim*dim-dim].blue + src[dim*dim-dim+1].blue) >> 2;

    // right_bottom
    dst[dim*dim-1].red = (src[dim*dim-dim-2].red + src[dim*dim-dim-1].red + src[dim*dim-2].red + src[dim*dim-1].red) >> 2;
    dst[dim*dim-1].green = (src[dim*dim-dim-2].green + src[dim*dim-dim-1].green + src[dim*dim-2].green + src[dim*dim-1].green) >> 2;
    dst[dim*dim-1].blue = (src[dim*dim-dim-2].blue + src[dim*dim-dim-1].blue + src[dim*dim-2].blue + src[dim*dim-1].blue) >> 2;

    // top
    sub_sum[0].red = (src[0].red + src[dim].red);
    sub_sum[0].green = (src[0].green + src[dim].green);
    sub_sum[0].blue = (src[0].blue + src[dim].blue);

    sub_sum[1].red = (src[1].red + src[dim+1].red);
    sub_sum[1].green = (src[1].green + src[dim+1].green);
    sub_sum[1].blue = (src[1].blue + src[dim+1].blue);

    sub_sum[2].red = (src[2].red + src[dim+2].red);
    sub_sum[2].green = (src[2].green + src[dim+2].green);
    sub_sum[2].blue = (src[2].blue + src[dim+2].blue);
    
    dst[1].red = (sub_sum[0].red+sub_sum[1].red+sub_sum[2].red) / 6;
    dst[1].green = (sub_sum[0].green+sub_sum[1].green+sub_sum[2].green) / 6;
    dst[1].blue = (sub_sum[0].blue+sub_sum[1].blue+sub_sum[2].blue) / 6;

    for(i=2; i<dim-1; ++i)
    {
        sub_sum[0]=sub_sum[1];
        sub_sum[1]=sub_sum[2];

        sub_sum[2].red = (src[i].red + src[dim+i].red);
        sub_sum[2].green = (src[i].green + src[dim+i].green);
        sub_sum[2].blue = (src[i].blue + src[dim+i].blue);

        pos = i-1;

        dst[pos].red = (sub_sum[0].red+sub_sum[1].red+sub_sum[2].red) / 6;
        dst[pos].green = (sub_sum[0].green+sub_sum[1].green+sub_sum[2].green) / 6;
        dst[pos].blue = (sub_sum[0].blue+sub_sum[1].blue+sub_sum[2].blue) / 6;
    }

    // bottom
    sub_sum[0].red = (src[dim*dim-dim-dim].red + src[dim*dim-dim].red);
    sub_sum[0].green = (src[dim*dim-dim-dim].green + src[dim*dim-dim].green);
    sub_sum[0].blue = (src[dim*dim-dim-dim].blue + src[dim*dim-dim].blue);

    sub_sum[1].red = (src[dim*dim-dim-dim+1].red + src[dim*dim-dim+1].red);
    sub_sum[1].green = (src[dim*dim-dim-dim+1].green + src[dim*dim-dim+1].green);
    sub_sum[1].blue = (src[dim*dim-dim-dim+1].blue + src[dim*dim-dim+1].blue);

    sub_sum[2].red = (src[dim*dim-dim-dim+2].red + src[dim*dim-dim+2].red);
    sub_sum[2].green = (src[dim*dim-dim-dim+2].green + src[dim*dim-dim+2].green);
    sub_sum[2].blue = (src[dim*dim-dim-dim+2].blue + src[dim*dim-dim+2].blue);
    
    dst[dim*dim-dim+1].red = (sub_sum[0].red+sub_sum[1].red+sub_sum[2].red) / 6;
    dst[dim*dim-dim+1].green = (sub_sum[0].green+sub_sum[1].green+sub_sum[2].green) / 6;
    dst[dim*dim-dim+1].blue = (sub_sum[0].blue+sub_sum[1].blue+sub_sum[2].blue) / 6;

    for(i=dim*dim-dim+2; i<dim*dim-1; i++)
    {
        ii = i+1;
        sub_sum[0]=sub_sum[1];
        sub_sum[1]=sub_sum[2];

        sub_sum[2].red = (src[ii-dim].red + src[ii].red);
        sub_sum[2].green = (src[ii-dim].green + src[ii].green);
        sub_sum[2].blue = (src[ii-dim].blue + src[ii].blue);

        dst[i].red = (sub_sum[0].red+sub_sum[1].red+sub_sum[2].red) / 6;
        dst[i].green = (sub_sum[0].green+sub_sum[1].green+sub_sum[2].green) / 6;
        dst[i].blue = (sub_sum[0].blue+sub_sum[1].blue+sub_sum[2].blue) / 6;
    }

    // left
    sub_sum[0].red = (src[0].red + src[1].red);
    sub_sum[0].green = (src[0].green + src[1].green);
    sub_sum[0].blue = (src[0].blue + src[1].blue);

    sub_sum[1].red = (src[dim].red + src[dim+1].red);
    sub_sum[1].green = (src[dim].green + src[dim+1].green);
    sub_sum[1].blue = (src[dim].blue + src[dim+1].blue);

    sub_sum[2].red = (src[dim+dim].red + src[dim+dim+1].red);
    sub_sum[2].green = (src[dim+dim].green + src[dim+dim+1].green);
    sub_sum[2].blue = (src[dim+dim].blue + src[dim+dim+1].blue);
    
    dst[1].red = (sub_sum[0].red+sub_sum[1].red+sub_sum[2].red) / 6;
    dst[1].green = (sub_sum[0].green+sub_sum[1].green+sub_sum[2].green) / 6;
    dst[1].blue = (sub_sum[0].blue+sub_sum[1].blue+sub_sum[2].blue) / 6;

    for(i=dim*3; i<dim*dim-1; i+=dim)
    {
        sub_sum[0]=sub_sum[1];
        sub_sum[1]=sub_sum[2];

        sub_sum[2].red = (src[i].red + src[i+1].red);
        sub_sum[2].green = (src[i].green + src[i+1].green);
        sub_sum[2].blue = (src[i].blue + src[i+1].blue);

        pos = i - dim;

        dst[pos].red = (sub_sum[0].red+sub_sum[1].red+sub_sum[2].red) / 6;
        dst[pos].green = (sub_sum[0].green+sub_sum[1].green+sub_sum[2].green) / 6;
        dst[pos].blue = (sub_sum[0].blue+sub_sum[1].blue+sub_sum[2].blue) / 6;
    }

    // right
    sub_sum[0].red = (src[dim-2].red + src[dim-1].red);
    sub_sum[0].green = (src[dim-2].green + src[dim-1].green);
    sub_sum[0].blue = (src[dim-2].blue + src[dim-1].blue);

    sub_sum[1].red = (src[dim+dim-2].red + src[dim+dim-1].red);
    sub_sum[1].green = (src[dim+dim-2].green + src[dim+dim-1].green);
    sub_sum[1].blue = (src[dim+dim-2].blue + src[dim+dim-1].blue);

    sub_sum[2].red = (src[dim*3-2].red + src[dim*3-1].red);
    sub_sum[2].green = (src[dim*3-2].green + src[dim*3-1].green);
    sub_sum[2].blue = (src[dim*3-2].blue + src[dim*3-1].blue);
    
    dst[1].red = (sub_sum[0].red+sub_sum[1].red+sub_sum[2].red) / 6;
    dst[1].green = (sub_sum[0].green+sub_sum[1].green+sub_sum[2].green) / 6;
    dst[1].blue = (sub_sum[0].blue+sub_sum[1].blue+sub_sum[2].blue) / 6;

    for(i=dim*4-1; i<dim*dim-1; i+=dim)
    {
        sub_sum[0]=sub_sum[1];
        sub_sum[1]=sub_sum[2];

        sub_sum[2].red = (src[i-1].red + src[i].red);
        sub_sum[2].green = (src[i-1].green + src[i].green);
        sub_sum[2].blue = (src[i-1].blue + src[i].blue);

        pos = i - dim;

        dst[pos].red = (sub_sum[0].red+sub_sum[1].red+sub_sum[2].red) / 6;
        dst[pos].green = (sub_sum[0].green+sub_sum[1].green+sub_sum[2].green) / 6;
        dst[pos].blue = (sub_sum[0].blue+sub_sum[1].blue+sub_sum[2].blue) / 6;
    }

    // other
    for(i=1; i<dim-1; i++)
    {
        j = 1;
        sub_sum[0].red = (src[(i-1)*dim].red+src[i*dim].red+src[(i+1)*dim].red);
        sub_sum[0].green = (src[(i-1)*dim].green+src[i*dim].green+src[(i+1)*dim].green);
        sub_sum[0].blue = (src[(i-1)*dim].blue+src[i*dim].blue+src[(i+1)*dim].blue);

        sub_sum[1].red = (src[(i-1)*dim+1].red+src[i*dim+1].red+src[(i+1)*dim+1].red);
        sub_sum[1].green = (src[(i-1)*dim+1].green+src[i*dim+1].green+src[(i+1)*dim+1].green);
        sub_sum[1].blue = (src[(i-1)*dim+1].blue+src[i*dim+1].blue+src[(i+1)*dim+1].blue);

        sub_sum[2].red = (src[(i-1)*dim+2].red+src[i*dim+2].red+src[(i+1)*dim+2].red);
        sub_sum[2].green = (src[(i-1)*dim+2].green+src[i*dim+2].green+src[(i+1)*dim+2].green);
        sub_sum[2].blue = (src[(i-1)*dim+2].blue+src[i*dim+2].blue+src[(i+1)*dim+2].blue);

        dst[i*dim+1].red = (sub_sum[0].red+sub_sum[1].red+sub_sum[2].red) / 9;
        dst[i*dim+1].green = (sub_sum[0].green+sub_sum[1].green+sub_sum[2].green) / 9;
        dst[i*dim+1].blue = (sub_sum[0].blue+sub_sum[1].blue+sub_sum[2].blue) / 9;

        for(j=2; j<dim-1; ++j)
        {
            sub_sum[0]=sub_sum[1];
            sub_sum[1]=sub_sum[2];

            sub_sum[2].red = (src[(i-1)*dim+j].red+src[i*dim+j].red+src[(i+1)*dim+j].red);
            sub_sum[2].red = (src[(i-1)*dim+j].red+src[i*dim+j].red+src[(i+1)*dim+j].red);
            sub_sum[2].red = (src[(i-1)*dim+j].red+src[i*dim+j].red+src[(i+1)*dim+j].red);

            pos = j-1;
            dst[i*dim+pos].red = (sub_sum[0].red+sub_sum[1].red+sub_sum[2].red) / 9;
            dst[i*dim+pos].green = (sub_sum[0].green+sub_sum[1].green+sub_sum[2].green) / 9;
            dst[i*dim+pos].blue = (sub_sum[0].blue+sub_sum[1].blue+sub_sum[2].blue) / 9;
        }
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
}

