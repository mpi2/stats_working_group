\name{PhenStat-NEWS}
\title{PhenStat News}
\section{CHANGES IN VERSION 2.3.1}{
  \subsection{NEW FUNCTIONALITY}{        
        \itemize{
            \item New output format is added: JSON format. In odred to get output in JSON 
            format apply JSONOutput with the PhenTestResult object as an argument.
        }
    }
    

  \subsection{COMPATIBILITY ISSUES}{
      
        \itemize{            
            \item The Box-Cox transformation is switched off by default. In order 
            to swithc in on the \code{testDataset()} function's argument “transformValues" 
            has to be set to TRUE.            
        }
    }  
}

\section{CHANGES IN VERSION 2.1.3}{
    
    \subsection{NEW FUNCTIONALITY}{        
        \itemize{
            \item Box-Cox transformation is applied when needed for the continues data. 
            Transformation values are shown in output functions. Genotype effects are converted back into 
            original scale if transformation was applied.
            \item New statistical method is implemented for the categorical data - Logistic Regression.
            \item New function is created \code{boxplotSexGenotypeBatchAdjusted()} that allows to evaluate Batch effect.
            \item Major architectural changes: classes PhenList and PhenTestResult are S4 classes now. 
            Different new methods are provided for the classes.
            \item PhenTestResult object contains dataset \code{analysedDataset()} that was used for the analysis including original values, transformed values (if transformation has been applied) and values adjusted for batch effect.
            
        }
    }
    \subsection{COMPATIBILITY ISSUES}{    
        \itemize{
            \item The \code{vectorOutput()} has new element for transformation values which have increased its length. 
            
        }
    }
}
\section{CHANGES IN VERSION 2.0.2}{
    
    \subsection{BUG FIXES}{
        
        \itemize{
            
            \item Column names in PhenList object.
             The check for the existing column names (Genotype, Batch, Sex, Weight) 
             is added before the renaming of the columns. 
            \item 'generateGraphs' function is back to work
            \item 'summaryOutput' function bug in FE framework is solved
            
        }
    }
}
\section{CHANGES IN VERSION 2.0.1}{
    
    \subsection{NEW FEATURES}{
        
        \itemize{
            
            \item Two new statistical methods implemented RR (Reference Range Plus) and TF (Time Fixed Effect) that can
            be called from \code{testDataset()} function with argument method set to RR and TF correspondingly.
            
            \item Additional measure of biological effect added: Percentage Change to \code{summaryOutput()} 
            and \code{vectorOutput} functions.
            
            \item New function is added to suggest analysis paths for dataset: \code{recommendMethod()}
            
            \item Function to implement the RR method: \code{RRTest()}.
            
            \item Functions to implement the TF method: \code{TFDataset()} creates dataset suitanle for TF, 
            \code{startTFModel()}, \code{finalTFModel} create model and fit it.
            
            \item Changes in function \code{summaryOutput()} - additional information and clearer layout.
            
        }
    }
    
    \subsection{COMPATIBILITY ISSUES}{
        
        \itemize{
            
            \item The \code{vectorOutput()} had additional elements which have increased its length. 
            \item The function \code{boxplotSexGenotypeBatch()} has been deprecated and replace with 
            \code{scatterplotSexGenotypeBatch()}. 
            \item Additional argument (phenotypeThreshold) with default value 0.01 has been added 
            to the \code{summaryOutput()} function
            
        }
    }
    
    \subsection{BUG FIXES}{
        
        \itemize{
            
            \item Number of critical bugs fixed. These are detailed at 
            \url{https://github.com/mpi2/stats_working_group/issues?q=is%3Aissue+is%3Aclosed} 
            GitHub software tracking tools.
            
        }
    }
  }
}