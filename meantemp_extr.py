#!/usr/bin/env python3
"""
Author : Emmanuel Gonzalez, Michele Cosi
Note   : Parts of this code was initially developed by the AgPipeline and TERRA-REF teams.
Date   : 2020-08-05
Purpose: Mean temp extraction using cv2 and stats
"""

import argparse
import cv2
import sys
import csv
import os
import uuid
import glob
import random
import statistics
import json
import pandas as pd
import tifffile as tifi
import matplotlib.pyplot as plt
import numpy as np
from osgeo import gdal
from PIL import Image
from scipy import stats
from scipy.signal import find_peaks

# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='extract mean temperature from plots (.tif); returns a .csv file per plot',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dir',
                        metavar='str',
                        help='Input directory of plot directories')

    parser.add_argument('-o',
                        '--outdir',
                        help='Output directory',
                        metavar='str',
                        type=str,
                        default='mean_temp_out')

    parser.add_argument('-O',
                        '--outname',
                        help='Output filename',
                        metavar='str',
                        type=str,
                        default='mean_temp')

    parser.add_argument('-g',
                        '--geo',
                        help='GeoJSON of plots',
                        type=str,
                        required=True)

    parser.add_argument('-d',
                        '--date',
                        help='processing date',
                        metavar='date of flir files you are processing',
                        type=str,
                        default='No date set')

    args = parser.parse_args()

    if '/' not in args.dir:
        args.dir = args.dir + '/'
    if '/' not in args.outdir:
        args.outdir = args.outdir + '/'

    return args

    # --------------------------------------------------
def get_trt_zones():
    trt_zone_1 = []
    trt_zone_2 = []
    trt_zone_3 = []

    for i in range(3, 19):
        for i2 in range(2, 48):
            plot = f'MAC_Field_Scanner_Season_10_Range_{i}_Column_{i2}'
            #print(plot)
            trt_zone_1.append(str(plot))

    for i in range(20, 36):
        for i2 in range(2, 48):
            plot = f'MAC_Field_Scanner_Season_10_Range_{i}_Column_{i2}'
            #print(plot)
            trt_zone_2.append(str(plot))

    for i in range(37, 53):
        for i2 in range(2, 48):
            plot = f'MAC_Field_Scanner_Season_10_Range_{i}_Column_{i2}'
            #print(plot)
            trt_zone_3.append(str(plot))

    return trt_zone_1, trt_zone_2, trt_zone_3

# --------------------------------------------------
def find_trt_zone(plot_name):
    trt_zone_1, trt_zone_2, trt_zone_3 = get_trt_zones()
    #print(trt_zone_1)

    if plot_name in trt_zone_1:
        trt = 'treatment 1'

    elif plot_name in trt_zone_2:
        trt = 'treatment 2'

    elif plot_name in trt_zone_3:
        trt = 'treatment 3'

    else:
        trt = 'border'

    return trt

# --------------------------------------------------
def get_genotype(plot, geo):
    with open(geo) as f:
        data = json.load(f)

    for feat in data['features']:
        if feat.get('properties')['ID']==plot:
            genotype = feat.get('properties').get('genotype')

    return genotype

# --------------------------------------------------
def main():

    args = get_args()
    temp_dict = {}
    temp_cnt = 0
    img_list = glob.glob(f'{args.dir}/*/*_ortho.tif', recursive=True)
    #print(img_list)

    for one_img in img_list:
        
        #### 1.1 PEAKS PROCESS
        temp_cnt += 1
        date = args.date
        plot_raw = one_img.split('/')[-2]
        genotype = get_genotype(plot_raw, args.geo)
        plot_name = '_'.join(plot_raw.split(' '))
        trt_zone = find_trt_zone(plot_name)
        #print(f'{plot_name}')
        #print(f'{trt_zone}\n')
        print(f'Processing {plot_raw}')

        g_img = gdal.Open(one_img)
        a_img = g_img.GetRasterBand(1).ReadAsArray()
        m = stats.mode(a_img)
        mode, count = m
        peak = mode[0][0:5].mean()
        temp = peak - 273.15

        a_img[a_img > peak] = np.nan
        mean_tc_peaks = np.nanmean(a_img) - 273.15

        #### 1.2 CV2 PROCESS

        gdal_in = gdal.Open(one_img)
        gdal_raster = gdal_in.GetRasterBand(1)
        im = gdal_raster.ReadAsArray()

        # if im is not None:
        normed = cv2.normalize(im, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)

        # Colour image gray
        im_color = cv2.applyColorMap(normed, cv2.COLORMAP_BONE)
        kernel = np.ones((1,1), np.uint8)

        #  Erode, dilate, and blur image for accurate contour detection
        img_erosion = cv2.erode(im_color, kernel, iterations=1)
        img_dilation = cv2.dilate(img_erosion, kernel, iterations=1)
        blur = cv2.GaussianBlur(img_dilation, (5,5), 1)           
        
        # Change colorspace and carry out edge detection for contour identification
        c_img = cv2.cvtColor(blur, cv2.COLOR_RGB2GRAY)
        kernel = np.ones((5,5), np.uint8)
        #edges = cv2.Canny(c_img,110,141)
        #res = cv2.morphologyEx(edges, cv2.MORPH_CLOSE, kernel)
        
        # Find max peak (odd) and constant
        #Const = 0
        peaks = cv2.minMaxLoc(c_img)
        if (peaks[1] % 2) == 0:
            o_peaks = ((peaks[1]-1)/2)+((peaks[1]+1)/2)+1
        else:
            o_peaks = ((peaks[1]-1)/2)+((peaks[1]+1)/2)
        if peaks[1] >= 250:
            Const = 62

        if  250 > peaks[1] >= 240:
            Const = 42

        if 240 > peaks[1] >= 230:
            Const = 32
        
        #maxpeak = peaks[1]
        # print(f'peaksloc: {peaks}\n') #minVal, maxVal, minLoc, maxLoc
        # print(f'globalmax: {peaks[1]}\n')
        # print(f'constant: {Const}\n')
        # print(f'oddity: {o_peaks}\n')
        try:

            gaus = cv2.adaptiveThreshold(c_img, int(peaks[1]), cv2.ADAPTIVE_THRESH_MEAN_C, 
                                            cv2.THRESH_BINARY, int(o_peaks),int(Const))
            
            kernel = np.ones((30, 30), np.uint8)
            closing = cv2.morphologyEx(gaus, cv2.MORPH_OPEN, kernel)
            
            contours, hierarchy = cv2.findContours(closing, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)           
            
            # Create mask and setup necessary variables for data collection
            canvas = im.copy()
            areas = {}
            num_cnt = 0
            shape = str(np.shape(canvas)).replace('(','').replace(')','')
            mask = np.zeros(im.shape, np.uint8)

            # Apply mask and extract name, meantemp, coords
            for c in contours:
                cnt = contours[0]
                M = cv2.moments(cnt)
                area = cv2.contourArea(c)

                if int(area) in range(10, 100000):
                    num_cnt += 1

                    # Create bounding box centered around contour
                    x,y,w,h = cv2.boundingRect(c)
                    cx = x+w//2
                    cy = y+h//2
                    center = (cx, cy)
                    dr = 30
                    cv2.rectangle(canvas, (int(cx-dr), int(cy-dr)), (int(cx+dr), int(cy+dr)), (0, 0, 255), 2)

                    # Draw the contours 
                    cv2.drawContours(im, contours, -1, (0,0,0), -1)

                    # Collect pixel coordinates and convert into lat, lon
                    pix_coord = str(center).strip('()')
                    x_pix, y_pix = pix_coord.split(',')
                    ds = gdal.Open(one_img)
                    c, a, b, f, d, e = ds.GetGeoTransform()
                    lon = a * int(x_pix) + b * int(y_pix) + a * 0.5 + b * 0.5 + c
                    lat = d * int(x_pix) + e * int(y_pix) + d * 0.5 + e * 0.5 + f
                    gps_coords = f'{lat}, {lon}'
                    
                    
                    # Crop each contour for individual plant temp measurements 
                    cropped = im[cy-dr:cy+dr, cx-dr:cx+dr]
                    imarray = np.array(cropped)
                    imarray[imarray == 0] = np.nan
                    #stdev = np.std(cropped)
                    mean_tc_cv2 = np.nanmean(imarray) - 273.15
                    print(f'Temp: {mean_tc_cv2}\n')
        except:
            pass               
    
        temp_dict[temp_cnt] = {
            'date': date,
            'treatment': trt_zone,
            'plot': plot_name,
            'genotype': genotype,
            'plot_temp': temp,
            'mean_temp_peaks': mean_tc_peaks,
            'mean_temp_cv2': mean_tc_cv2
            }
    df = pd.DataFrame.from_dict(temp_dict, orient='index', columns=['date',
                                                                    'treatment',
                                                                    'plot',
                                                                    'genotype',
                                                                    'plot_temp',
                                                                    'mean_temp_peaks',
                                                                    'mean_temp_cv2'])
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    df.to_csv(os.path.join(args.outdir, args.outname + '.csv'), index=False)   

    print(f'Done. Check outputs in {args.outdir}')

# --------------------------------------------------
if __name__ == '__main__':
    main()
