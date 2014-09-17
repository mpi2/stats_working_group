\name{PhenStat-NEWS}
\title{PhenStat News}

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
    
    \subsection{COMPTABILITY ISSUES}{
        
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