library(stringr)

extract.deepest.taxonomy <- function(fullTax) {
    taxonomy = ""
    
    if(str_detect(fullTax, "\\.s__")) {
        tmp = str_split(fullTax, "\\.s__")[[1]]
        taxonomy = tmp[length(tmp)]
    } else {
        if(str_detect(fullTax, "\\.g__")) {
            tmp = str_split(fullTax, "\\.g__")[[1]]
            taxonomy = tmp[length(tmp)]
        } else {
            if (str_detect(fullTax, "\\.f__")) {
                tmp = str_split(fullTax, "\\.f__")[[1]]
                taxonomy = tmp[length(tmp)]
            } else {
                if (str_detect(fullTax, "\\.o__")) {
                    tmp = str_split(fullTax, "\\.o__")[[1]]
                    taxonomy = tmp[length(tmp)]
                } else {
                    if (str_detect(fullTax, "\\.c__")) {
                        tmp = str_split(fullTax, "\\.c__")[[1]]
                        taxonomy = tmp[length(tmp)]
                    } else {
                        if (str_detect(fullTax, "\\.p__")) {
                            tmp = str_split(fullTax, "\\.p__")[[1]]
                            taxonomy = tmp[length(tmp)]
                        } else {
                            if (str_detect(fullTax, "d__")) {
                                tmp = str_split(fullTax, "d__")[[1]]
                                taxonomy = tmp[length(tmp)]
                            } else {
                                taxonomy = fullTax
                            }      
                        }
                        
                    }
                }
                
            }
        }
    }
    taxonomy = str_replace_all(taxonomy, "\\.__", "")
    return(taxonomy)
}