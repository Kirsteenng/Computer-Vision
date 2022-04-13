#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    // TODO Fill this in
    //checking boundaries
    x = x >= im.w ? im.w-1: x ;
    x = x < 0 ? 0 : x;
    y =  y >= im.h ? im.h-1: y ;
    y =  y  <0 ?0: y ;
    return im.data[x + (y*im.w) + (c*im.w*im.h)];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // TODO Fill this in
    // check coordinates has gone beyond boundaries
    if (x < im.w && y < im.h && x >=0 && y>=0){
    im.data[x + (y*im.w) + (c*im.w*im.h)] = v;  
    }
    

}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    // TODO Fill this in
    memcpy(copy.data, im.data, im.w*im.h*im.c*sizeof(im.w));
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
    float r,g,b;
    for(int width = 0; width <= im.w -1; width++){
        for (int height = 0; height<= im.h -1; height++){
                r = get_pixel(im,width,height,0);
                g = get_pixel(im,width,height,1);
                b = get_pixel(im,width,height,2);
                gray.data[width + height*im.w] = 0.299*r + 0.587*g + 0.114*b;
            }
        }

    return (gray);
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in
    for (int width = 0; width <= im.w; width++){
        for (int height = 0; height <= im.h ; height++){
            set_pixel(im, width, height, c, (get_pixel(im,width,height,c) + v));
        }
    }
}

void clamp_image(image im)
{
    // TODO Fill this in
    float pixel;
    for (int width =0; width < im.w; width ++){
        for(int height = 0; height < im.h; height ++){
            for(int channel = 0; channel < im.c;channel++){
                pixel = get_pixel(im,width,height,channel);
                pixel = pixel < 0 ? 0: pixel;
                pixel = pixel > 1 ? 1: pixel;
                set_pixel(im,width,height,channel,pixel);

            }

        }
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    // TODO Fill this in
    float value, hue, hue_p, saturation,r,g,b,min_value,max_value,c;
    int width,height;

    for (width = 0; width < im.w ; width++){
        for (height = 0; height < im.h; height ++){
            r = get_pixel(im,width,height,0);
            g = get_pixel(im,width,height,1);
            b = get_pixel(im,width,height,2);

            min_value = three_way_min(r,g,b);
            max_value = three_way_max(r,g,b);

            value = max_value;
            c = value - min_value;

            if (r == 0 && g == 0 && b ==0){
                saturation = 0;
            }
            else{
                 saturation = c/value;
            }

            if ( c == 0)
            {
                hue_p = 0;
            }
            else if(value == r){
                hue_p = (g - b)/c;
            }
            
             else if(value == g){
                hue_p = ((b - r)/c) + 2;
            }
             else if(value == b){
                hue_p = ((r - g)/c) + 4;
            }

            if(hue_p < 0){
                hue = (hue_p/6) + 1;
            }
            else{
                hue = hue_p /6;
            }

            // setting R channel with H, the G channel with S and the B channel with V.
            set_pixel(im,width,height,0,hue);
            set_pixel(im,width,height,1,saturation);
            set_pixel(im,width,height,2,value);
            }
            
        }

    }

void hsv_to_rgb(image im)
{
    // TODO Fill this in
    float h, F, P, Q,T,H,S,V;
    int h_i;
    int width, height;
    for (width = 0; width < im.w ; width++){
        for (height = 0; height < im.h; height ++){
            H = get_pixel(im,width,height,0);
            S = get_pixel(im,width,height,1);
            V = get_pixel(im,width,height,2);
            h = H * 6;
            h_i = floor(h);
            F = h - h_i;
            P = V * (1-S);
            Q = V * (1 - F*S);
            T = V * (1 - (1 - F) * S);

            switch(h_i){
                case 0:
                    set_pixel(im,width,height,0,V);
                    set_pixel(im,width,height,1,T);
                    set_pixel(im,width,height,2,P);
                    break;

                case 1:
                    set_pixel(im,width,height,0,Q);
                    set_pixel(im,width,height,1,V);
                    set_pixel(im,width,height,2,P);
                    break;

                case 2:
                    set_pixel(im,width,height,0,P);
                    set_pixel(im,width,height,1,V);
                    set_pixel(im,width,height,2,T);
                    break;

                case 3:
                    set_pixel(im,width,height,0,P);
                    set_pixel(im,width,height,1,Q);
                    set_pixel(im,width,height,2,V);
                    break;

                case 4:
                    set_pixel(im,width,height,0,T);
                    set_pixel(im,width,height,1,P);
                    set_pixel(im,width,height,2,V);
                    break;
                
                case 5:
                   set_pixel(im,width,height,0,V);
                    set_pixel(im,width,height,1,P);
                    set_pixel(im,width,height,2,Q);
                    break;


            }

        }

    }
}
// Extra credit
void scale_image(image im, int c, float v){
    int width, height;
    float target_pixel;
    for (width = 0; width <= im.w; width ++){
        for (height = 0; height <= im.h; height ++){
            set_pixel(im, width, height,c, (get_pixel(im,width,height,c) * v));
        }
    }
}
