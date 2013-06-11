outputLength <- function(result)
# Technical results needed for the output
# The output lengths in the vector form depends of effects that are included into the model
{
    if(!is(result,"PhenTestResult")) {
        stop("Please provide result paramtere as PhenTestResult object")
    }
    
    numberofgenders=result$numberGenders
    keep_weight <- result$weightEffect
    keep_gender <- result$genderEffect
    keep_interaction <- result$interactionEffect
    keep_batch <- result$batchEffect
    keep_equalvar <- result$varianceEffect
    equation <- result$equation
    
    table_length <- NA
    
    if (equation=="withWeight"){   
        if(numberofgenders==2){
            if((keep_gender && keep_interaction)|(!keep_gender && keep_interaction)){
                table_length=5
            }else if(keep_gender && !keep_interaction){
                table_length=4
            }else {
                table_length=3
            }  
        }else{
            table_length=3
        }
    }
    # Eq.1
    else{
        if(numberofgenders==2){
            if((keep_gender && keep_interaction)|(!keep_gender && keep_interaction)){
                table_length=4
            }else if(!keep_gender && !keep_interaction){
                table_length=2
            }else{
                table_length=3
            }  
        }else{
            table_length=2
        }
    } 
    
    return (table_length)
}


