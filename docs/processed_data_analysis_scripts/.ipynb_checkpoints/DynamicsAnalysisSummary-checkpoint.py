# library packages
import os
import pandas as pd
import numpy as np


def fileLoad(folder_input, file_input):
    file_list = []
    for file in file_input:
        file_path = os.path.join(folder_input, file)
        try:
            df = pd.read_csv(file_path, index_col=None, encoding='utf-8')
        except:
            df =pd.read_csv(file_path, index_col=None, encoding='cp1252')
        file_list.append(df)
    return file_list 

parameter_files = (fileLoad(folder, fileTable))

def cellVolumeLength():
    cell_volume = parameter_files[0]
    cell_vol = cell_volume['Vol (µm³)']
    cell_volume_std = cell_volume['std_Vol']
    
    # calculate for cell length
    cell_length = cell_vol**(1/3)
    cell_volume['Cell length (µm)'] = cell_length
    
    # calculate for cell length propagation error
    Error_prop = abs((1/3) * cell_vol**(-2/3) * cell_volume_std)
    cell_volume['Error propagation'] = Error_prop
    
    return cell_volume

def chromosomeDistance(chromosome_file_number):
    chromosome_to_chromosome_distance = parameter_files[chromosome_file_number]
    chromosome_mean = pd.DataFrame()
    chromosome_std = pd.DataFrame()
    for col_chromosome in chromosome_to_chromosome_distance:
        '''calculate for mean chromosome-to-chromosome distance for each embryonic cell stage'''
        chrom_temp_val_mean = pd.DataFrame.from_dict({f"{col_chromosome}": np.mean(chromosome_to_chromosome_distance[col_chromosome])},
                                                     orient='index',
                                                     columns=['Final chromosome-to-chromosome (µm)'])
        
        '''calculate for standard deviation chromosome-to-chromosome distance for each embryonic cell stage'''
        chrom_temp_val_std = pd.DataFrame.from_dict({f"{col_chromosome}": np.std(chromosome_to_chromosome_distance[col_chromosome])},
                                                     orient='index',
                                                     columns=['Final chromosome-to-chromosome std'])
        
        # new chromosome-to-chromosome mean and standard deviation tables
        chromosome_mean = pd.concat([chromosome_mean, chrom_temp_val_mean], axis=0)
        chromosome_std = pd.concat([chromosome_std, chrom_temp_val_std], axis=0)
        
    # new table
    chromosome_distance = pd.concat([chromosome_mean, chromosome_std], axis=1)
    return chromosome_distance

def poleDistance():
    pole_to_pole_distance = parameter_files[2]
    poles_mean = pd.DataFrame()
    poles_std = pd.DataFrame()
    for col_pole in pole_to_pole_distance:
        '''calculate for mean pole-to-pole distance for each embryonic cell stage'''
        pole_temp_val_mean = pd.DataFrame.from_dict({f"{col_pole}": np.mean(pole_to_pole_distance[col_pole])}, orient='index', columns=['Final pole-to-pole (µm)'])
        
        '''calculate for standard deviation of the pole-to-pole distance for each embryonic cell stage'''
        pole_temp_val_std = pd.DataFrame.from_dict({f"{col_pole}": np.std(pole_to_pole_distance[col_pole])}, orient='index', columns=['Final pole-to-pole std'])
        
        # new pole-to-pole mean and standard deviation tables
        poles_mean = pd.concat([poles_mean, pole_temp_val_mean], axis=0)
        poles_std = pd.concat([poles_std, pole_temp_val_std], axis=0)
        
    # new table 
    pole_distance = pd.concat([poles_mean, poles_std], axis=1)
    return pole_distance

def poleChromosomeDistance():
    pole_to_chromosome_distance = parameter_files[3]
    pole_chrom_mean = pd.DataFrame()
    pole_chrom_std = pd.DataFrame()
    for col_pole_chrom in pole_to_chromosome_distance:
        '''calculate for mean pole-to-chromosome distance for each embryonic cell stage'''
        pole_chrom_temp_mean = pd.DataFrame.from_dict({f"{col_pole_chrom}": np.mean(pole_to_chromosome_distance[col_pole_chrom])},
                                                      orient='index', columns=['Final pole-to-chromosome (µm)'])
        
        '''calculate for standard pole-to-chromosome distance for each embryonic cell stage'''
        pole_chrom_temp_std = pd.DataFrame.from_dict({f"{col_pole_chrom}": np.std(pole_to_chromosome_distance[col_pole_chrom])},
                                                      orient='index', columns=['Final pole-to-chromosome std'])
        
        # new pole-to-chromosome mean and standard deviation tables
        pole_chrom_mean = pd.concat([pole_chrom_mean, pole_chrom_temp_mean], axis=0)
        pole_chrom_std = pd.concat([pole_chrom_std, pole_chrom_temp_std], axis=0)
        
    # new tables
    polesChroms_distance = pd.concat([pole_chrom_mean, pole_chrom_std], axis=1)
    return polesChroms_distance

def metaphasePoleDistance():
    metaphase_pole_to_pole_distance = parameter_files[4]
    metaphase_pole_mean = pd.DataFrame()
    metaphase_pole_std = pd.DataFrame()
    for col_metaphase in metaphase_pole_to_pole_distance:
        '''calculate for mean metaphase pole-to-pole distance for each embryonic cell stage'''
        meta_pole_temp_mean = pd.DataFrame.from_dict({f"{col_metaphase}": np.mean(metaphase_pole_to_pole_distance[col_metaphase])},
                                                      orient='index', columns=['Metaphase pole-to-pole (µm)'])
        
        '''calculate for standard deviation metaphase pole-to-pole distance for each embryonic cell stage'''
        meta_pole_temp_std = pd.DataFrame.from_dict({f"{col_metaphase}": np.std(metaphase_pole_to_pole_distance[col_metaphase])},
                                                      orient='index', columns=['Metaphase pole-to-pole std'])
        
        # new metaphase pole-to-pole mean and standard deviation tables
        metaphase_pole_mean = pd.concat([metaphase_pole_mean, meta_pole_temp_mean], axis=0)
        metaphase_pole_std = pd.concat([metaphase_pole_std, meta_pole_temp_std], axis=0)
        
    # new tables
    metaphase_distance = pd.concat([metaphase_pole_mean, metaphase_pole_std], axis=1)
    return metaphase_distance

def elongationRatePoles():
    elongation_speed_pole_to_pole = parameter_files[5]
    elongation_speed_mean = pd.DataFrame()
    elongation_speed_std = pd.DataFrame()
    for col_elongation in elongation_speed_pole_to_pole:
        '''calculate for mean pole-to-pole distance elongation rate for each embryonic cell stage'''
        elong_speed_temp_mean = pd.DataFrame.from_dict({f"{col_elongation}": np.mean(elongation_speed_pole_to_pole[col_elongation])},
                                                      orient='index', columns=['Elongation rate pole-to-pole (µm/min)'])
        
        '''calculate for standard deviation pole-to-pole distance elongation rate for each embryonic cell stage'''
        elong_speed_temp_std = pd.DataFrame.from_dict({f"{col_elongation}": np.std(elongation_speed_pole_to_pole[col_elongation])},
                                                      orient='index', columns=['Elongation rate pole-to-pole std'])
        
        # new pole-to-pole elongation rate mean and standard deviation tables
        elongation_speed_mean = pd.concat([elongation_speed_mean, elong_speed_temp_mean], axis=0)
        elongation_speed_std = pd.concat([elongation_speed_std, elong_speed_temp_std], axis=0)
        
    # new tables
    elongation_speed = pd.concat([elongation_speed_mean, elongation_speed_std], axis=1)
    return elongation_speed

def segregationSpeedChromosomes():
    segregation_speed_chromosome_to_chromosome = parameter_files[6]
    segregation_speed_mean = pd.DataFrame()
    segregation_speed_std = pd.DataFrame()
    for col_segregation in segregation_speed_chromosome_to_chromosome:
        '''calculate for mean chromosome-to-chromosome distance segregation speed for each embryonic cell stage'''
        elong_speed_temp_mean = pd.DataFrame.from_dict({f"{col_segregation}": np.mean(segregation_speed_chromosome_to_chromosome[col_segregation])},
                                                      orient='index', columns=['segregation speed chromosome-to-chromosome (µm/min)'])
        
        '''calculate for standard deviation chromosome-to-chromosome distance segregation speed for each embryonic cell stage'''
        elong_speed_temp_std = pd.DataFrame.from_dict({f"{col_segregation}": np.std(segregation_speed_chromosome_to_chromosome[col_segregation])},
                                                      orient='index', columns=['segregation speed chromosome-to-chromosome std'])
        
        # new chromosome-to-chromosome segregation speed mean and standard deviation tables
        segregation_speed_mean = pd.concat([segregation_speed_mean, elong_speed_temp_mean], axis=0)
        segregation_speed_std = pd.concat([segregation_speed_std, elong_speed_temp_std], axis=0)
        
    # new tables
    segregation_speed = pd.concat([segregation_speed_mean, segregation_speed_std], axis=1)
    return segregation_speed

def initialPoleDistance():
    initial_pole_distance = parameter_files[7]
    initial_mean = pd.DataFrame()
    initial_std = pd.DataFrame()
    for col_initial in initial_pole_distance:
        '''calculate for mean initial pole-to-pole distance for each embryonic cell stage'''
        initial_temp_val_mean = pd.DataFrame.from_dict({f"{col_initial}": np.mean(initial_pole_distance[col_initial])}, orient='index', columns=['Initial pole-to-pole (µm)'])
        
        '''calculate for std initial pole-to-pole distance for each embryonic cell stage'''
        initial_temp_val_std = pd.DataFrame.from_dict({f"{col_initial}": np.std(initial_pole_distance[col_initial])}, orient='index', columns=['Initial pole-to-pole std'])
        
        # new initial pole-to-pole distance mean and standard deviation tables
        initial_mean = pd.concat([initial_mean, initial_temp_val_mean], axis=0)
        initial_std = pd.concat([initial_std, initial_temp_val_std], axis=0)
        
    # new tables
    initial_poles = pd.concat([initial_mean, initial_std], axis=1)
    return initial_poles

def spindleDynamicsValues():
    cell_volume_length = cellVolumeLength(cell_volume_file_number)
    final_chromosome_length = chromosomeDistance()
    final_poles_length = poleDistance()
    final_pole_chromosome = poleChromosomeDistance()
    metaphase_length = metaphasePoleDistance()
    elongation = elongationRatePoles()
    segregation = segregationSpeedChromosomes()
    initial_pole_length = initialPoleDistance()
    
    # table with spindle dynamics values
    newTable = pd.concat([initial_pole_length, final_poles_length, final_chromosome_length, final_pole_chromosome, metaphase_length, elongation, segregation], axis=1).reset_index(drop=True)
    new_df = pd.concat([cell_volume_length, newTable], axis=1)
    
    # save table as csv file
    csv_new_df = new_df.to_csv(os.path.join(save_files, 'Dynamics_analysis_summary ().csv'), index=False, encoding='cp1252')
    
    return csv_new_df