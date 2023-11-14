# packages 
import os
import glob
import matplotlib
import numpy as np
import pandas as pd
from skimage import io
from pathlib import Path
from scipy.signal import find_peaks
from matplotlib import pyplot as plt
from scipy.signal import argrelextrema
from skimage.io import imread, imshow, imsave
from skimage.transform import rescale, resize,downscale_local_mean

######################################################################################################################

def KyMeas(folder, ext):
    """
    take a parameter folder and return the list of filenames within the folder 
    which has the ending ext
    """
    Exp = []
    
    # iterate through the files in the folder
    for file in os.listdir(folder):
        _, file_ext = os.path.splitext(file) #get the file extension
        
        # get the path of each file within the folder and append to the list 
        if ext in file_ext:
            file_path = os.path.join(folder, file)
            file_path = file_path.replace("\\", "/")
            Exp.append(file_path)
            
    return Exp 

def find_maxima(img1, prominence_threshold, channel=1): 
    
    """detecting the 2 maximum peaks for every timepoint in each of the images"""
    maxima = []
    
    for j in range(img1.shape[1]):
        intensity = img1[channel, j]
        x = intensity
        """prominence measures the peak that stands out from the surrounding baseline of signal intensity"""
        peaks, _ = find_peaks(x, prominence = (prominence_threshold, None)) #   
        #plt.plot([j]*len(peaks), peaks, ".")   
        maxima.append(peaks)
    
    return maxima

def filterMaxima(list_of_maxima_pairs, distance_threshold):
    """
    this function account for 2 maxima for every time point (iteration)
    * the distance threshold is the maximum distance travelled between 2 timepoints. 
    """
    
    # filter out the third maxima 
    temp = list_of_maxima_pairs[0]
    temp.sort() # rearrange maxima in each frame in ascending order by sort function
    filtered_maxima = [temp[-2:]] # take the last 2 maxima in each frame
    #filtered_maxima = [temp[:2]]
    #filtered_maxima = [temp[1:]]
    
    for n in range(1, len(list_of_maxima_pairs)):
        first_maxima_set = list_of_maxima_pairs[n]
        second_maxima_set = filtered_maxima[n-1]
        
        # create a list of new array in the filtered array
        new_maxima = []
        for first_maximum in first_maxima_set:
            #print(first_maximum) """check"""
            for second_maximum in second_maxima_set:
                distance = abs(first_maximum - second_maximum)
               
                if distance < distance_threshold:
                    new_maxima.append(first_maximum)
                    break
        
        new_maxima.sort()
        filtered_maxima.append(new_maxima[-2:])
        #filtered_maxima.append(new_maxima[0:2])
        #filtered_maxima.append(new_maxima[1:])
    
    return filtered_maxima

def findChromosomeMaxima(maximas_centrosome, maximas_chromosome):
        '''
        This function is used to find the chromosome maxima.
        The idea is to find the maximas within the boundary of the centrosome maxima. 
        '''
        maxChromosome_between_maxCentrosome = []
        for maxima_centrosome, maxima_chromosome in zip(maximas_centrosome, maximas_chromosome):
        
            start_centrosome = min(maxima_centrosome)
            end_centrosome = max(maxima_centrosome)
            maxima_chromosome = np.asarray(maxima_chromosome)
            mask = np.logical_and(maxima_chromosome > start_centrosome, maxima_chromosome < end_centrosome)

            maxChromosome_between_maxCentrosome.append(maxima_chromosome[mask])
        return maxChromosome_between_maxCentrosome

def plotMaxima(list_of_maxima_pair):
    plt.figure()
    for timepoint, peaks in enumerate(list_of_maxima_pair):
        plt.plot([timepoint]*len(peaks), peaks, ".")
    plt.show()    
    
def kymoDistance(new_maxima):
    frames = []
    distances = []
    pixel_size = 0.1
    for frame, i in enumerate(new_maxima):
        if len(i) == 2:            
            a = i[0]
            b = i[1]
            
            # absolute distance between 2 points
            distance = round((abs(a-b))*pixel_size, 3)
            distances.append(distance)
            frames.append(frame+1)
    
    return frames, distances

def plotDistances(filename, folder, frames, distances):
    plt.figure()
    a = frames
    b = distances
    plt.xlabel("time (sec)")
    plt.ylabel("distance (μm)")
    plt.plot(a, b)
    plt.title(filename)
    filename = os.path.join(folder, filename.split("/")[-1]).replace("\\", "/")
    plt.savefig(filename)
    plt.show()

def savecsvfile(filename, folder, frames, distances):
    pole_to_pole_distance = {"frame": frames, "distance": distances}
    df = pd.DataFrame(pole_to_pole_distance)
    df_files = os.path.join(folder, filename.split("/")[-1]).replace("\\", "/")
    df.to_csv(df_files, index=False)
    
def guessThreshold(input_image):
    
    """calculate the average threshold of all the maxima and output the 2 maxima with the average threshold"""
    
    count = 0
    total_threshold = 0
    for threshold in range(0, 10000, 10):
        maxima_x = find_maxima(input_image, threshold)
        maxima_x = filterMaxima(maxima_x, 50)
        
        count_maxima_first_frame = len(maxima_x[0])
        count_maxima_sixth_frame = len(maxima_x[5])
        #print(count_maxima_first_frame, count_maxima_sixth_frame, threshold, "count_maxima")
        
        if count_maxima_first_frame == 2:
            total_threshold = total_threshold + threshold
            count = count + 1 
    
    Avg_threshold = total_threshold/count
    print("Average_threshold:", Avg_threshold)
        
    return Avg_threshold

def ChromosomeMaximaQuality(list_maxima, maximum_frame = 20):
        count_correct_detection = 0
        for num, maxima in enumerate(list_maxima):
            if num >= maximum_frame:
                break
            #print(num, maxima)
            
            number_maxima = len(maxima)
            if number_maxima >0 and number_maxima < 3:
                #print(number_maxima)
                count_correct_detection=count_correct_detection +1
            
        return count_correct_detection/maximum_frame

    
"""" Reading out the positions of the centrosomes and chromosome maxima and saving as a dataframe"""    
def reformatPositionsPoles(maxima_positions):
    minima_position = [np.min(f) for f in maxima_positions]
    maxima_position = [np.max(f) for f in maxima_positions]
    position_table = pd.DataFrame({
        'pole_1':minima_position,
        'pole_2':maxima_position,
    })
    return position_table

def reformatPositionsChromosome(maxima_positions):
    minima_position = [np.min(f) for f in maxima_positions]
    maxima_position = [np.max(f) for f in maxima_positions]
    position_table = pd.DataFrame({
        'chromosome_1':minima_position,
        'chromosome_2':maxima_position,
    })
    return position_table