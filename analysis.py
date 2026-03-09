#importing important packages
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import Counter
from matplotlib.patches import ConnectionPatch

#class for visualisation/summary of the publications in my systematic review
class PublicationsSummary: 
    '''Class specific to the publications csv,'''

    def __init__(self,publications_df,geo_df):
        self.publications=publications_df
        self.geo=geo_df
        self.pub=self.publications.copy()
        self.pub.loc[self.pub['Animal model']!='human',"Animal model"]='Other'
  

    def human_vs_nonhuman(self):
        '''Function which counts how many clock publications were human specific and how many 
        were specific to other animals. 
        Function returns a pie chart'''

        if 'Animal model' not in self.publications.columns:
            raise ValueError("Missing 'Animal model' column.")

        #counting the amount of human and non-human clocks 
        count_human=self.pub['Animal model'].value_counts().get('human')
        count_other=self.pub['Animal model'].value_counts().get('Other')
        return(count_human,count_other)

        #Naming variables for the pie chart
        data=[float(count_human),float(count_other)]
        labels=['Human','Other']

        #Creating the pie chart
        plt.figure(figsize=(10,6))
        #defining Seaborn platte color to use 
        palette_color=sns.color_palette('dark')
        #plotting data on chart
        plt.pie(data,labels=labels,autopct='%.0f%%',startangle=90,textprops={'fontsize':14})
        plt.show()

    def all_species(self):
        '''Function which returns a dictionary of counts for each species'''
        list_species=list(self.publications['Animal model'])
        animals_count={} #final (unsorted) where the count for each species will be saved
        sorted_animals_count={} #final sorted dictionary

        #refining the list
        final_species_list=[]
        for animal in list_species:
            animal_str=str(animal)
            final_species_list.extend(animal_str.split(','))

        #Creating a dictionary which counts how many times each animal model appears 
        for animal in final_species_list:
            animals_count[animal]=final_species_list.count(animal)
        
        #Sorting the dictionary 
        sorted_animals_count=dict(sorted(animals_count.items(),key=lambda x:x[1],reverse=True))

        return(sorted_animals_count)
    
    def epi_clocks_over_time(self):
        '''Function which returns a histogram of the amount of epigenetic clock publications over the years '''
        #Specifiying multiple figures in one
        plt.rcParams.update({'font.size': 20})
        fig, ax=plt.subplots()
        fig.set_size_inches(25,15)
        #first histogram
        ax=sns.histplot(data=self.pub,x='Year',hue='Animal model', multiple='stack')
        ax.set_ylabel('Count of Epigenetic Clocks')
        #specify that the line plot will have the same x axis as the histplot
        ax2=ax.twinx()
        sns.lineplot(x='Year',y='Total GEO Datasets',data=self.geo,ax=ax2,color='red',label='GEO datasets')
        ax2.legend(loc='upper left')
        plt.savefig('epi_clocks_time.pdf')
        plt.show()
        






class MethylationMeasurement:
    '''Specific to epigenetic studies, for methylation measurement methods only'''

    def __init__(self,publications_df):
        self.clocks=publications_df
        self.human=self.clocks[self.clocks['Animal model']=='human']
        self.nonhuman=self.clocks[~self.clocks['Animal model'].str.contains('human',case=False,na=False)]


    def cleaning(self,method):
        '''Function which removes unecessary punctuations and replaces terms with other equivalent terms'''
        if str(method).lower()=='nan':
            return None
        method = str(method)
        method = method.replace('.', '')
        method = method.replace('?', '')
        method = method.replace(' assay', '')
        

        #Replacement of equivalent terms
        replacements = {
            'Illumina Infinium MethylationEPIC BeadChip': 'array',
            '450K': 'array',
            '27K': 'array',
            'Beadchip': 'array',
            '850K': 'array',
            'EPIC': 'array',
            'Illumina HumanMethylation arrays': 'array',
            'mouse methylation beadchip':'array',
            'HorvathMammalMethyl40': 'array',
            'HorvathMammalMethylChip40': 'array',
            'Pyrosequencing': 'pyrosequencing',
            'Qiagen PyroMark':'pyrosequencing',
            'pyrosequencings':'pyrosequencing',
            'bisulfite polymerase chain reaction (PCR) sequencing': 'bisulfite PCR',
            'Targeted bisulphite sequencing': 'Targeted bisulfite sequencing',
            'Targeted bisulfite Sequencing': 'Targeted bisulfite sequencing',
            'scRRBS': 'RRBS',
            'scWGBS': 'WGBS'
        }

        for old,new in replacements.items():
            method=method.replace(old,new)
        
        return method
    
    def methods_dictionary(self,methods_list):
        '''Function which outputs a dictionary of methylation measurement counts'''
        split_methods=[x.strip() for item in methods_list for x in item.split(',')]
        cleaned_methods=[self.cleaning(method) for method in split_methods if method]
        method_dic=dict(Counter(cleaned_methods))

        return method_dic




    
    def human_methylation_methods(self):
        '''Function that will call all the previous functions to output a pie chart of the different methylation measurement methods'''
        methylation_method=list(self.human['Methylation Measurement Method'])
        method_dic=self.methods_dictionary(methylation_method)
        #need to re-arrange the dictionary
        dic_for_pie={'EpiTYPER': 0, 'pyrosequencing': 0, 'array': 0, 'other': 0}
        for key, value in method_dic.items():
            if key in ('EpiTYPER', 'pyrosequencing', 'array'):
                dic_for_pie[key] += value
            else:
                dic_for_pie['other'] += value
        
        
        #here we need to correct the method_dic because it has included too many arrays. In some clocks they've used datasets coming from multiple arrays. 
        keywords=['450K','27K','BeadChip','Beadchip','850K','EPIC','Illumina HumanMethylation arrays','Illumina Infinium MethylationEPIC BeadChip']
        count=0
        for tech in methylation_method:
            if any(word in tech for word in keywords):
                count=count+1
        corrected_method_dic={'EpiTYPER': dic_for_pie['EpiTYPER'], 'pyrosequencing': dic_for_pie['pyrosequencing'], 'array': count, 'other': dic_for_pie['other']}
        #pie=self.piechart(method_dic=corrected_method_dic,array_dict=array_dict_human,explode_list=[0,0,0.1,0],array_index=2,angle_constant=180)
        return corrected_method_dic,method_dic


    
    def nonhuman_methylation_methods(self):
        methylation_method=list(self.nonhuman['Methylation Measurement Method'])
        method_dic=self.methods_dictionary(methylation_method)
        final_dic={'Other':0}
        for key,value in method_dic.items():
            if key in ('Oxford Nanopore Technologies','Bisulfite-converted restriction site-associated DNA sequencing','multiplex PCR','bisulfite PCR sequencing','Oxford Nanopore Technology','bisulfite amplicon high-throughput sequencing'):
                final_dic['Other']+=value
            else:
                final_dic[key]=value
        final_dic.pop('illumina miseq',None)
        return(final_dic)
        
       












        
        
        
        
        
