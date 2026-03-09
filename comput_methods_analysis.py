#### Package for the analysis of computational methods used in epigenetic clocks such as Imputation, feature selection and age prediction




#importing important packages
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict,Counter

class ImputationMethods: ##Issue here as 'no, removed missing is considered an imputation method, although it isn't'
    '''Class which gives all the relevant information about imputation 
    across the different clocks reviewed'''
    
    def __init__(self,epi_clocks):
        self.clocks=epi_clocks
        #cleaning steps
        #Drop the unecessary columns in epi_clocks
        #self.clocks.drop(columns=['train Spearman Coefficient','train MAD','train MAE','train Pearson','test Spearman Coefficient','test MAD','all Spearman Coefficient','MAD all','test Pearson correlation','Test R2','MAE TEST','MedAE TEST','Median error test','Test MAD2','MSE','RMSE'],inplace=True)
        #replacing average with mean because it means the same thing
        self.clocks.loc[self.clocks['Imputation']=='average','Imputation']='mean'
        self.clocks.loc[self.clocks['Imputation']=='Y','Imputation']=None
        self.clocks.loc[self.clocks['Imputation']=='Multivariate Imputation by Chained Equations','Imputation']='mice'
        self.clocks.loc[self.clocks['Feature selection method']=='none','Feature selection method']=None
        
        #differentiate human and non human clocks
        self.human=self.clocks[self.clocks['Animal model']=='human']
        self.nonhuman=self.clocks[~self.clocks['Animal model'].str.contains('human',case=False,na=False)]
        self.imput=self.clocks[~self.clocks['Imputation'].isna()]
        #Dictionary of dataframes and their names
        self.dataframes={'all the studies':self.clocks, 
                         'human studies':self.human,
                         'non-human animal studies':self.nonhuman}
        
        
        self.assays={'pyrosequencing':['pyrosequencing','Pyrosequencing','pyromark'],
               'sequencing':['RRBS','WGBS','BBA-seq','targeted bisulfite sequencing','sequencing','seq','miseq'],
               'illumina array':['EPIC','450K','850K','27K','Illumina HumanMethylation arrays',
               'Illumina Infinium MethylationEPIC BeadChip','HorvathMammalMethylChip40','Beadchip'],
               'epityper':['EpiTYPER','epityper']
               }
        
        #assays are more detailed here 
        self.detailed_assays = {'pyrosequencing':['pyrosequencing','Pyrosequencing','pyromark'],
                'targeted sequencing':['BBA-seq','targeted bisulfite sequencing','miseq','time-seq','amplicon','Bisulfite-converted restriction site-associated DNA sequencing'],
               'genome-wide sequencing':['RRBS','WGBS','sequencing'],
               'illumina array':['EPIC','450K','850K','27K','Illumina HumanMethylation arrays',
               'Illumina Infinium MethylationEPIC BeadChip','HorvathMammalMethylChip40','Beadchip'],
               'epityper':['EpiTYPER','epityper']
               }
        
    def imputation_present_absent(self):
        '''Function which creates barchart illustrating the proportion of studies that specified using imputat
        imputation method and those that didn't''' #this f
        for name,dataframe in self.dataframes.items():
            imputations=list(dataframe['Imputation'].notnull().value_counts(True))
            return dataframe['Imputation'].notnull().value_counts()
            labels=['No imputation method','Imputation method']
            plt.bar(labels,imputations,color=['steelblue','lightgray'])
            plt.title(f'Proportion of Imputation methods in studies for {name}')
            plt.show()
            
            
    
    
 


class FeatureSelectionMethods:
    '''Class to analyze the different feature selection methods used in the epigenetic clocks'''
    
    def __init__(self,epi_clocks:pd.DataFrame):
        self.clocks=epi_clocks
        self.clocks.loc[self.clocks['Feature selection method']=='none','Feature selection method']=None
        #differentiate human and non human clocks
        self.human=self.clocks[self.clocks['Animal model']=='human']
        self.nonhuman=self.clocks[~self.clocks['Animal model'].str.contains('human',case=False,na=False)]
        #list of the feature selection methods used in human clocks
        self.fslist_human=list(self.human['Feature selection method'])
        #creating categories for feature selection methods of human clocks 
        self.fsmethods_human = {
    "Pre-selection Filters": [
        "psp", "psr", "cpgs with reliable ddpcr assays", "coverage 50x", "iqr>0.1",
        "average of single-cpg age estimates","filtering"
    ],
    "Correlation-Based Methods": [
        "pearson correlation", "spearman correlation", "spearman rank correlation", 
        "correlation with age", "correlation with logarithmic age", "pearson log age", 
        "spearman absolute correlation", "weighted correlation network analysis (wgcna)", 
        "q-values", "absolute correlation 0.2"
    ],
    "Other Linear Methods": [
        "linear regression", "stepwise regression", "forward stepwise regression", 
        "stepwise forward feature selection", "stepwise multiple linear regression",
        "multivariate regression analysis", "linear coefficient analysis", 
        "linear mixed effect model (ewas)", "multivariable stepwise linear regression with bic",
    ],
    "Linear Penalized Regression (lasso or elastic net)": [
        "elastic net", "lasso", "lasso regression", "causality-enriched elastic net model",
        "elastic net bootstrap", "stability selection framework(with lasso)"
    ],
    "Tree-Based": [
        "random forest regression", "lightgbm"
    ],
    "Neural Network": [
        "deep feature selection gradient based saliency",
    ],
    "Other ML / Algorithmic Methods": [
        "recursive feature elimination", "%-rfe to 100", "%-rfe to 1500", "%-rfe to 10000",
        "kbest25", "kbest2000", "boruta", "sfm(extra trees)", "sfm(elas)", "genetic algorithm",
        "mutual information", "sure independence screening", "importance ranking", 
        "frequency prioritization", "fdr<0.05", "fdr", "pavlidis template matching",'super learner','pca'
    ],
    "Causal Inference": [
        "epigenome-wide mendelian randomization", "causality-enriched elastic net model"
    ],
    "None / Not Specified": [
        "none"
    ]
}

        
        #list of the feature selection methods used in nonhuman clocks
        self.fslist_nonhuman=list(self.nonhuman['Feature selection method'])
        #creating categories for feature selection methods of nonhuman clocks
        self.fsmethods_nonhuman=self.fsmethods_nonhuman = {
    "Pre-selection Filters": [
        "hm", "psr"
    ],
    "Correlation-Based Methods": [
        "pearson correlation", "pearson correlation coefficient", "spearman correlation", 
        "correlation analysis", "backup site", "linear regression, pearson correlation",
        "p-value filtering", "median z statistic of a correlation test"
    ],
    "Other Linear  Methods": [
        "linear regression", "stepwise linear regression",'bayesian generalized linear models'
    ],
    "Linear Penalized Regression (elastic net or lasso)": [
        "elastic net", "lasso"
    ],

    
    "Other ML Methods": [
        "methylation relationship matrice 1", "methylation relationship matrice 2", 
        "unsupervised clustering", "ranked intersection algorithm"
    ]
}

        
        
    # def categorizing_feature_selection(self,feature_selection_list,fsmethods:dict):
    #     '''Function which takes in a list of feature selection method and a dictionary of feature selection categories.
    #     It will count how many times each category is present in the list of feature selection methods and will output that as a dictionary of counts'''
    #     #defining variable
    #     categorized_feature_selection={}
    #     for method in feature_selection_list:
    #         if not method is None and isinstance(method,str): #if method is none, handling missing data here
    #             methods=[m.strip().lower() for m in method.replace('\n',',').split(',')] #removing trailing white spaces and put everything lower case 
    #             for m in methods:
    #                 found=False #Need to add reset for each method, this specifies that it hasn't found the category yet. or else without the found it will stop at the first category
    #                 for category,keywords in fsmethods.items():
    #                     if any(keyword.lower() in m for keyword in keywords):
    #                         categorized_feature_selection[category]=categorized_feature_selection.get(category,0)+1 #dictionary empty at the start 
    #                         found=True
    #     #sorting everything in reverse order for graph 
    #     #sorted_feature_selection_counts={k:v for k, v in sorted(categorized_feature_selection.items(), key=lambda item:item[1],reverse=True)}
    #     sorted_feature_selection_counts=feature_selection_list

    #     return sorted_feature_selection_counts
    
    def categorizing_feature_selection(self, feature_selection_list, fsmethods: dict):
        categorized_feature_selection = {}

        for method in feature_selection_list:
            if method is None or not isinstance(method, str) or not method.strip():
                continue

            tokens = [m.strip().lower()
                    for m in method.replace('\n', ',').split(',')
                    if m.strip()]

            #categories hit by THIS clock
            categories_for_this_clock = set()

            for tok in tokens:
                for category, keywords in fsmethods.items():
                    if any(keyword.lower() in tok for keyword in keywords):
                        categories_for_this_clock.add(category)

            #increment each category ONCE per clock
            for cat in categories_for_this_clock:
                categorized_feature_selection[cat] = categorized_feature_selection.get(cat, 0) + 1

        #sort high->low (optional)
        sorted_feature_selection_counts=dict(sorted(categorized_feature_selection.items(), key=lambda x: x[1], reverse=True))
        return sorted_feature_selection_counts

    

    
    def plotting_fs_fequencies(self, sorted_feature_selection_counts: dict,title:str):
        '''Function which plots a bar chart of the different feature selection categories with count labels'''
        plt.figure(figsize=(12, 6))
        categories = list(sorted_feature_selection_counts.keys())
        counts = list(sorted_feature_selection_counts.values())

        plt.barh(categories, counts)

        # Add count labels next to bars
        for i, count in enumerate(counts):
            plt.text(count + 0.5, i, str(count), va='center')  # Adjust 0.5 if needed for spacing

        # Labels and title
        plt.xlabel("Count")
        plt.ylabel("Feature Selection Categories")
        plt.title("Frequency of Feature Selection Categories")

        # Show the plot
        plt.gca().invert_yaxis()  # Highest count at the top
        plt.tight_layout()
        plt.savefig(f'feature_selection_{title}.pdf')
        plt.show()
 
        
    def human_feature_selection(self):
        sorted_feature_selection_counts=self.categorizing_feature_selection(self.fslist_human,self.fsmethods_human)
        #return self.human['Feature selection method'].str.contains('spearman',case=False).sum()
        #return sorted_feature_selection_counts
        plot=self.plotting_fs_fequencies(sorted_feature_selection_counts=sorted_feature_selection_counts,title='human')
        return plot
    
    def nonhuman_feature_selection(self):
        sorted_feature_selection_counts=self.categorizing_feature_selection(self.fslist_nonhuman,self.fsmethods_nonhuman)
        #return self.nonhuman['Feature selection method'].str.contains('spearman',case=False).sum()
        #return sorted_feature_selection_counts
        plot=self.plotting_fs_fequencies(sorted_feature_selection_counts=sorted_feature_selection_counts,title='non_human')
        return plot
    

class AgePredictionMethods:
    '''Class to analyze the different feature selection methods used in the epigenetic clocks'''
    def __init__(self,epi_clocks):
        self.clocks=epi_clocks
        self.human=self.clocks[self.clocks['Animal model']=='human']
        self.nonhuman=self.clocks[~self.clocks['Animal model'].str.contains('human',case=False,na=False)]
        
        #list of human models
        self.humanmodels=list(self.human['Model'])
        #list of human catgeories
        self.humancategories = {
    "Linear Penalized Regression (lasso, ridge or elastic net)": [
        "elastic net", "lasso", "ridge regression", "cox with elastic net"
    ],
    "Other Linear Models": [
        "linear regression", "multiple linear regression", "multivariate linear regression",
        "predefined linear regression", "linear-regression",
        "stepwise regression", "stepwise multiple regression", "multivariable linear model",
        "multivariable stepwise linear regression", "multivariable model", "multivariate regression model",'simulated annealing approach was used to train linear models optimizing for model accuracy, significance of the correlation between delta age and lifestyle/health factors, and model complexity'
    ],
    "Quantile Regression": [
        "quantile regression", "multivariate quantile regression"
    ],
    "Support Vector Machines": [
        "support vector machine", "support vector machine (radial kernel)",
        "support vector regression", "svr (polynomial)", "svmr"
    ],
    "Neural Networks / Deep Learning": [
        "neural network model", "deep neural network", "artificial neural network (multilayer perceptron)",
        "back propagation neural network", "generalized regression neural network", "resnet",
        "perceptron prediction model based on the channel attention mechanism"
       
    ],
    "Tree-Based Models": [
        "random forest regression", "xgboost", "lightgbm", "catboost",
        "gradient boosting regressor model", "gbr"
    ],
    
    "Probabilistic / Likelihood-Based Models": [
        "gaussian process regression",
        "maximum likelihood estimation +LOWESS using a binomial distribution to model the probability of observing the counts of methylated and unmethylated reads at each CpG site"
    ],
    "Other / Hybrid ": [
        "mmda", "mqr", "wenda", "quantile regression support vector machine", "quantile regression neural network", "generalized additive model", "super learner","quadratic regression model"
    ]
}

        
        #list of non-human models
        self.nonhumanmodels = list(self.nonhuman['Model'])

        #list of nonhuman categories
        self.nonhumancategories = {
    "Linear Penalized Regression (lasso, ridge or elastic net)": [
        "elastic net", "penalized regression", "elastic net + age transformation",
        "lasso", "lasso regression", "ridge regression",
        "multinomial logistic regression with elastic net"
    ],
    "Other Linear Models": [
        "single regression", "multiple linear regression", "multivariate regression",
        "multivariate linear regression", "reverse regression"
    ],
    "Support Vector Machines": [
        "support vector"
    ],
    "Tree-Based Models": [
        "random forest"
    ],
    "Probabilistic / Likelihood-Based Models": [
        "maximum likelihood model"
    ],
    "Other / Mixed Models": [
        "meblup", "epigenetic pacemaker", "average of single-cpg age estimates", "random forest,elastic net"
    ]
}
        
    def categorizing_models(self,raw_models:list,categories:dict):
        '''Function which takes in one list of models and one dictionary which links models to lists.
        Function will categorize each model present in the clocks and will output a dictionary of counts.
        '''
        #assign categories to models
        category_counts = defaultdict(int)
        uncategorized = []
        for model in raw_models:
            model_clean = model.lower().strip()
            found = False
            for cat, keywords in categories.items():
                if any(keyword.lower() in model_clean for keyword in keywords):
                    category_counts[cat] += 1
                    found = True
                    break
            if not found:
                uncategorized.append(model)

        #sorting for plotting
        sorted_counts = dict(sorted(category_counts.items(), key=lambda item: item[1], reverse=True))
        return sorted_counts
        
    def plotting(self,sorted_counts:dict,title:str):
        '''Function which takes in list and outputs a bar plot of model categories'''
        
        plt.figure(figsize=(10, 6))
        plt.barh(list(sorted_counts.keys()), list(sorted_counts.values()), color="steelblue")
        
                # Add count labels next to bars
        for i, count in enumerate(sorted_counts.values()):
            plt.text(count + 0.5, i, str(count), va='center')
        plt.xlabel("Count")
        plt.title("Model Types Used in Epigenetic Clocks")
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(f'age_prediction_{title}.pdf')
        plt.show()
        
    def human_models(self):
        '''Function for human clocks'''
        sorted_counts=self.categorizing_models(raw_models=self.humanmodels,categories=self.humancategories)
        #return sorted_counts
        title='human'
        plotting=self.plotting(sorted_counts=sorted_counts,title=title)
        
        return plotting
    
    def nonhuman_models(self):
        '''Function for nonhuman clocks'''
        sorted_counts=self.categorizing_models(raw_models=self.nonhumanmodels,categories=self.nonhumancategories)
        #return sorted_counts
        title='non_human'
        plotting=self.plotting(sorted_counts=sorted_counts,title=title)
        return plotting
        
        
        
    
        
    
#####################################################################################################
############################################# 3 Functions below are only used within other functions ####################################          
             
        
#    def count_methods(self,df):
#         '''Function which counts the prevalence of each methylation measurement method'''
#         counts= {key:0 for key in self.assays}
#         for assay in df['Methylation Measurement Method']:
#             assay_lower=assay.lower()
#             for method,keywords in self.assays.items():
#                 if any(keyword.lower() in assay_lower for keyword in keywords):
#                     counts[method]+=1
#                     break
#         return counts   
       

#     def platform_type(self,assay,assay_dict):
#             '''Function used in the function imputation_per_methylation_measurement.
#             It categorizes each methylation measurement method into types'''
#             assay_lower=assay.lower()
#             for type,keywords in assay_dict.items():
#                 if any(keyword.lower() in assay_lower for keyword in keywords):
#                     return(type)
#             return('other')
        
#     def imputation_type(self,imputation):
#         imput_type={'Mean':['mean'],
#             'KNN':['knn'],
#             'Other simple':['filled with 0','median'],
#             'Other ML':['mice','imputPCA()','linear regression','incorporated in model','softImpute ALS option','elastic net']
#             }
#         imputation_lower=imputation.lower()
#         for type,keywords in imput_type.items():
#             if any(keyword.lower() in imputation_lower for keyword in keywords):
#                 return(type)
               
# #########################################################################################################
          
#     def methylation_measurement_imputation_or_not(self):
#         counts_epi=self.count_methods(df=self.clocks)
#         counts_imput=self.count_methods(df=self.imput)
#         #percentages for the graph
#         percentages={key:counts_imput[key]/counts_epi[key]*100 for key in counts_epi.keys()}
#         #Graph of methods 
#         methods=list(counts_imput.keys())
#         count_total_values=list(counts_epi.values())
#         count_imput_values=list(counts_imput.values())
#         count_not_imput=[counts_epi[method]-counts_imput[method] for method in methods]
#         #bar width for the graph
#         plt.figure(figsize=(10,6))
#         bar_width=0.4

#         #Percentages
#         p_imput=[counts_imput[method]/counts_epi[method]*100 for method in methods]
#         p_not=[(counts_epi[method]-counts_imput[method])/counts_epi[method]*100 for method in methods]
#         percentages=[p_imput,p_not]


#         #creating an array to create numerical positions for the x-axis categories in the bar chart
#         x=np.arange(len(methods))
#         graph1=plt.bar(x-bar_width/2,count_imput_values,bar_width,label='Imputation Used',alpha=0.7,)
#         graph2=plt.bar(x+bar_width/2,count_not_imput,bar_width,label="No Imputation Used",alpha=0.7)
#         graphs=[graph1,graph2]

#         plt.xlabel('Methylation Measurement Method')
#         plt.ylabel('Number of Epigenetic Clocks')
#         plt.title('Use of Imputation by Methylation Measurement Method')
#         plt.xticks(ticks=x, labels=methods, rotation=45)


#         z=np.arange(2)

#         # for i,(graph,p) in enumerate(zip(graphs,percentages)):
#         #     for j,bar in enumerate(graph):
#         #         width=bar_width
#         #         height=bar.get_height()
#         #         plt.text(bar.get_x()+width/2,height+1, f"{p[j]:.1f}%", ha='center', fontsize=10)
        
#         count_lists = [count_imput_values, count_not_imput]

#         for i, (graph, counts) in enumerate(zip(graphs, count_lists)):
#             for j, bar in enumerate(graph):
#                 width = bar_width
#                 height = bar.get_height()
#                 plt.text(bar.get_x() + width / 2, height + 1, f"{int(counts[j])}", ha='center', fontsize=10)


#         plt.legend()
#         plt.show()


#     def imputation_per_methylation_measurement(self):
#         '''Function which outputs a bar/histplot of the different imputation methods used for each methylation measurement methods.
#         Mainly a plotting function'''
#         self.imput['platform_type']=self.imput['Methylation Measurement Method'].apply(lambda x: self.platform_type(x,self.detailed_assays))
#         self.imput['imputation_type']=self.imput['Imputation'].apply(lambda x: self.imputation_type(imputation=x))
#         #Code for histplot
#         plt.figure(figsize=(16, 10))
#         sns.histplot(data=self.imput,x='platform_type',hue='imputation_type',multiple='stack')
#         plt.xlabel('Methylation Platform Type')
#         plt.xticks(rotation=75)
#         plt.show()
        