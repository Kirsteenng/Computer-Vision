#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

const int HIGHPASS[3][3] = {{0, -1, 0},
                            {-1, 4, -1},
                            {0, -1, 0}};

const int SHARPEN[3][3] = {{0, -1, 0},
                           {-1, 5, -1},
                           {0, -1, 0}};

const int EMBOSS[3][3] = {{-2, -1, 0},
                          {-1, 1, 1},
                          {0, 1, 2}};

/******************************** Resizing *****************************
  To resize we'll need some interpolation methods and a function to create
  a new image and fill it in with our interpolation methods.
************************************************************************/

float nn_interpolate(image im, float x, float y, int c)
{
    // Question: what is floating column and how does it relate to offset? 
    /***********************************************************************
      This function performs nearest-neighbor interpolation on image "im"
      given a floating column value "x", row value "y" and integer channel "c",
      and returns the interpolated value.
    ************************************************************************/
    
    return get_pixel(im, (int)round(x),(int)round(y),c);
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix the return line)
    /***********************************************************************
      This function uses nearest-neighbor interpolation on image "im" to a new
      image of size "w x h"
    ************************************************************************/
    float new_x, new_y;
    int width, height, depth;
    image result = make_image(w,h,im.c);
    float scale_x = (float) im.w/(float)w;
    float scale_y = (float) im.h/ (float)h;
    float offset = 0.5;
    for (width = 0; width < w; width ++){
      for(height = 0; height < h; height ++){
        for(depth = 0; depth < im.c; depth ++){
          new_x = (width + offset)*scale_x - offset;
          new_y = (height + offset)*scale_y - offset;
          set_pixel(result,width, height,depth,nn_interpolate(im, new_x, new_y,depth));
        }

      }
    }
    return result;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs bilinear interpolation on image "im" given
      a floating column value "x", row value "y" and integer channel "c".
      It interpolates and returns the interpolated value.
    ************************************************************************/
// x and y will be the calculated new coordinates

float p00,p01,p10,p11,weightage;
float Ap00,Ap01,Ap10,Ap11;

int floor_x = (int)floor(x);
int floor_y = (int)floor(y);
//printf("Cood (0,0) =   (%d, %d) \n",floor_x,floor_y);

float distx0 = x - floor_x;
float disty0 = y - floor_y;


p00 = get_pixel(im,floor_x,floor_y,c);
p01 = get_pixel(im,floor_x, floor_y + 1, c);
p10 = get_pixel(im,floor_x+1, floor_y, c);
p11 = get_pixel(im,floor_x+1, floor_y + 1, c);

Ap00 = distx0 * disty0;
Ap10 = (1-distx0) * disty0;
Ap01 = distx0 * (1-disty0);
Ap11 = (1-distx0) * (1-disty0);

weightage = (p00 * Ap11) + (p01* Ap10) + (p10 *Ap01) + (p11 * Ap00); 

return weightage;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    /***********************************************************************
      This function uses bilinear interpolation on image "im" to a new image
      of size "w x h". Algorithm is same as nearest-neighbor interpolation.
    ************************************************************************/
    image result = make_image(w,h,im.c);

    float new_x, new_y;
    int width, height, depth;
    
    float scale_x = (float) im.w/w;
    float scale_y = (float) im.h/h;
    float offset = 0.5;
    for (width = 0; width < w; width ++){
      for(height = 0; height < h; height ++){
        for(depth = 0; depth < im.c; depth ++){
          new_x = (width + offset)*scale_x - offset;
          new_y = (height + offset)*scale_y - offset;
          set_pixel(result,width, height,depth,bilinear_interpolate(im, new_x, new_y,depth));
        }
      }
    }

    return result;
}


/********************** Filtering: Box filter ***************************
  We want to create a box filter. We will only use square box filters.
************************************************************************/

void l1_normalize(image im)
{
    // TODO
    /***********************************************************************
      This function divides each value in image "im" by the sum of all the
      values in the image and modifies the image in place.
    ************************************************************************/
  
  int width, height, depth;
  float normalized_pixel, sum = 0;

   for (width = 0; width < im.w; width ++){
      for(height = 0; height < im.h; height ++){
        for(depth = 0; depth < im.c; depth ++){
        sum += get_pixel(im,width, height, depth);
        }
      }
   }
   for (width = 0; width < im.w; width ++){
      for(height = 0; height < im.h; height ++){
        for(depth = 0; depth < im.c; depth ++){
        normalized_pixel = get_pixel(im,width, height, depth)/sum;
        set_pixel(im, width,height,depth,normalized_pixel);
        }
      }
   }
}

image make_box_filter(int w)
{
    /***********************************************************************
      This function makes a square filter of size "w x w". Make an image of
      width = height = w and number of channels = 1, with all entries equal
      to 1. Then use "l1_normalize" to normalize your filter.
    ************************************************************************/
  image box_filter = make_image(w, w, 1);
  for (int k = 0; k < box_filter.c; k++) {
    for (int j = 0; j < box_filter.h; j++) {
      for (int i = 0; i < box_filter.w; i++) {
        set_pixel(box_filter, i, j, k, 1);
      }
    }
  }
  l1_normalize(box_filter);
  return box_filter;
}

float apply_filter(image im, image filter, int w, int h, int imC, int filterC) {
  float res = 0;
  for (int j = 0; j < filter.h; j++) {
    for (int i = 0; i < filter.w; i++) {
      int newH = h + j - filter.h / 2;
      int newW = w + i - filter.w / 2;
      res += get_pixel(im, newW, newH, imC) * get_pixel(filter, i, j, filterC);
    }
  }
  return res;
}
// preserve the channel
void preserve_filter_convolve(image im, image filter, image res, int filterC) {
  for (int k = 0; k < im.c; k++) {
    for (int j = 0; j < im.h; j++) {
      for (int i = 0; i < im.w; i++) {
        set_pixel(res, i, j, k, apply_filter(im, filter, i, j, k, filterC == -1 ? k : filterC));
      }
    }
  }
}
// channel = 1
void no_preserve_filter_convolve(image im, image filter, image res, int filterC) {
  for (int i = 0; i < im.w; i++) {
    for (int j = 0; j < im.h; j++) {
      float sum = 0;
      for (int k = 0; k < im.c; k++) {
        sum += apply_filter(im, filter, i, j, k, filterC == -1 ? k : filterC);
      }
      set_pixel(res, i, j, 0, sum);
    }
  }
}

image convolve_image(image im, image filter, int preserve)
{
   assert(im.c==filter.c || filter.c==1);
   float val=0;
   float filter_p, im_p;
   image new_im;
   int currw,currh;
   if (im.c==filter.c && preserve == 0)
   {
     printf("Entering case 1, filter.c = %d, im.c = %d, preserve = %d\n", filter.c, im.c,preserve);
    new_im = make_image(im.w, im.h, 1);
     for (int i=0;i<im.w;i++)
    {
        for(int j=0;j<im.h;j++)
        {
          val=0;
          for(int k=0;k<im.c;k++)
          {
            for (int l=0;l<filter.h;l++)
            {
              for (int w=0;w<filter.w;w++)
              {
                  filter_p = get_pixel(filter,w,l,k);
                  currw = (int) i-filter.w/2 + w;
                  currh = (int) j-filter.h/2 +l;
                  im_p = get_pixel(im,currw,currh,k);
                  val += filter_p * im_p;
              }          
            }
          }
          set_pixel(new_im,i,j,0, val);
        }
    }
   }
   else if (im.c==filter.c && preserve == 1)
   {
     printf("Entering case 2, filter.c = %d, im.c = %d, preserve = %d\n", filter.c, im.c,preserve);
     new_im = make_image(im.w, im.h, im.c);
     for (int k=0;k<im.c;k++)
      {
        for(int j=0;j<im.h;j++)
        {
          for(int i=0;i<im.w;i++)
          {
            val=0;
            for (int l=0;l<filter.h;l++)
            {
              for (int w=0;w<filter.w;w++)
              {
                  filter_p = get_pixel(filter,w,l,k);
                  currw = (int) i-filter.w/2 + w;
                  currh = (int) j-filter.h/2 +l;
                  im_p = get_pixel(im,currw,currh,k);
                  val += filter_p * im_p;
              }
            }
            set_pixel(new_im,i,j,k,(float)val);
          }
        }
        
      }
   }
  else if (im.c!=filter.c && preserve == 0)
   {//im.c != filter.c && preserve == 0,put eveything into one channel
   printf("Entering case 3, filter.c = %d, im.c = %d, preserve = %d\n", filter.c, im.c,preserve);

       new_im = make_image(im.w, im.h, 1);
        for (int i=0;i<im.w;i++)
        {
            for(int j=0;j<im.h;j++)
            {
              val=0;
              for(int k=0;k<im.c;k++)
              {
                for (int l=0;l<filter.h;l++)
                {
                  for (int w=0;w<filter.w;w++)
                  {
                  filter_p = get_pixel(filter,w,l,0);
                  currw = (int) i-filter.w/2 + w;
                  currh = (int) j-filter.h/2 +l;
                  im_p = get_pixel(im,currw,currh,k);
                  val += filter_p * im_p;
                  }
                }
              }
              set_pixel(new_im,i,j,0,val);
            }
          }
   }
     else
    {
      //im.c != filter.c && preserve == 1
    // printf("Entering case 4, filter.c = %d, im.c = %d, preserve = %d\n", filter.c, im.c,preserve);
     new_im = make_image(im.w, im.h, im.c);
      for (int k=0;k<im.c;k++)
        {
          for(int j=0;j<im.h;j++)
          {
            for(int i=0;i<im.w;i++)
            {
              val=0;
              for (int l=0;l<filter.h; l++)
              {
                for (int w=0;w<filter.w; w++)
                {
                  filter_p = get_pixel(filter,w,l,0);
                  currw = (int) i- ceil(filter.w/2) + w;
                  currh = (int) j-ceil(filter.h/2) +l;
                  im_p = get_pixel(im,currw,currh,k);
                  val += filter_p * im_p;
                }       
              }
              set_pixel(new_im,i,j,k,val);
            }
          }
        }  
    }
    return new_im;
}


image make_highpass_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with highpass filter values using image.data[]
    ************************************************************************/
    image result = make_image(3, 3, 1);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
      set_pixel(result, i, j, 0, HIGHPASS[i][j]);
      }
  }
  return result;
}

image make_sharpen_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with sharpen filter values using image.data[]
    ************************************************************************/
image result = make_image(3,3,1);
for (int i =0; i < 3; i++){
  for(int j = 0; j< 3; j++){
    set_pixel(result, i, j, 0, SHARPEN[i][j]);
  }
}
    return result;
}

image make_emboss_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with emboss filter values using image.data[]
    ************************************************************************/
    image result = make_image(3,3,1);
for (int i =0; i < 3; i++){
  for(int j = 0; j< 3; j++){
    set_pixel(result, i, j, 0, EMBOSS[i][j]);
  }
}
    return result;
}

// Question 2.3.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: 
//Preserve == 1; EMBOSS and SHARPEN. Because in both cases we are trying to enchance the edges of the image, the color should be preserved
//Preserve == 0; HIGHPASS. We are extracting the edge from the image using this filter so color detail does not matter. 
//

// Question 2.3.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: 
// Post processing: filters that use preserve because we might want to apply smoothening on the edges to make it more defined.
// We do not want to post process filters such as highpass that highlights with the edges because that would hinder us from identifying the edges.
//

float calc_pdf(float sigma, int x, int y){
  return  exp((-(x * x + y * y) / (2 * sigma * sigma))) / (TWOPI * sigma * sigma) ;
}
image make_gaussian_filter(float sigma)
{
    // TODO
    /***********************************************************************
      sigma: a float number for the Gaussian.
      Create a Gaussian filter with the given sigma. Note that the kernel size 
      is the next highest odd integer from 6 x sigma. Return the Gaussian filter.
    ************************************************************************/
    int dimension = (int) ceil(sigma) * 6;
    dimension  = (dimension % 2 == 0 ? dimension + 1 : dimension);

    image gaussian = make_image(dimension,dimension,1);
    for (int width = 0 ; width < gaussian.w; width ++){
      for (int height = 0 ; height < gaussian.h; height ++){
        set_pixel( gaussian, width , height, 0, calc_pdf(sigma,width - gaussian.w/2,height - gaussian.h/2)); // why need to divide by 2?
      }
    }
    l1_normalize(gaussian);// why need to normalise?
    return gaussian;
}

image add_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input images a and image b have the same height, width, and channels.
      Sum the given two images and return the result, which should also have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
  assert( a.w == b.w && a.h == b.h && a.c == b.c);
  image result = make_image(a.w, a.h,a.c);
  float pix_a, pix_b;
  int width, height, channel;
  for (channel = 0; channel < a.c ; channel++){
    for(height = 0; height < a.h ; height++){
      for(width = 0; width < a.w ; width++){
        pix_a = get_pixel(a, width, height, channel);
        pix_b = get_pixel(b, width, height, channel);
        set_pixel(result, width, height,channel, pix_a + pix_b);
      }
    }
  }
  return result;
}

image sub_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input image a and image b have the same height, width, and channels.
      Subtract the given two images and return the result, which should have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    return make_image(1,1,1);
}

image make_gx_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gx filter and return it
    ************************************************************************/
    return make_image(1,1,1);
}

image make_gy_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gy filter and return it
    ************************************************************************/
    return make_image(1,1,1);
}

void feature_normalize(image im)
{
    // TODO
    /***********************************************************************
      Calculate minimum and maximum pixel values. Normalize the image by
      subtracting the minimum and dividing by the max-min difference.
    ************************************************************************/
}

image *sobel_image(image im)
{
    // TODO
    /***********************************************************************
      Apply Sobel filter to the input image "im", get the magnitude as sobelimg[0]
      and gradient as sobelimg[1], and return the result.
    ************************************************************************/
    image *sobelimg = calloc(2, sizeof(image));
    return sobelimg;
}

image colorize_sobel(image im)
{
  // TODO
  /***********************************************************************
    Create a colorized version of the edges in image "im" using the 
    algorithm described in the README.
  ************************************************************************/
  return make_image(1,1,1);
}

// EXTRA CREDIT: Median filter

/*
image apply_median_filter(image im, int kernel_size)
{
  return make_image(1,1,1);
}
*/

// SUPER EXTRA CREDIT: Bilateral filter

/*
image apply_bilateral_filter(image im, float sigma1, float sigma2)
{
  return make_image(1,1,1);
}
*/
