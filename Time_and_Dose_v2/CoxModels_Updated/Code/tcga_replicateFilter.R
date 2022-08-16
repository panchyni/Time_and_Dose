#---------------------------------------------------
# Filter TCGA Replicate Samples
# Author: ShixiangWang <w_shixiang@163.com>
# ooooooo
# The filter rule following broad institute says:
#
# In many instances there is more than one aliquot for a given combination of individual, platform, and data type. However, only one aliquot may be ingested into Firehose. Therefore, a set of precedence rules are applied to select the most scientifically advantageous one among them. Two filters are applied to achieve this aim: an Analyte Replicate Filter and a Sort Replicate Filter.
# 
# Analyte Replicate Filter
# The following precedence rules are applied when the aliquots have differing analytes. For RNA aliquots, T analytes are dropped in preference to H and R analytes, since T is the inferior extraction protocol. If H and R are encountered, H is the chosen analyte. This is somewhat arbitrary and subject to change, since it is not clear at present whether H or R is the better protocol. If there are multiple aliquots associated with the chosen RNA analyte, the aliquot with the later plate number is chosen. For DNA aliquots, D analytes (native DNA) are preferred over G, W, or X (whole-genome amplified) analytes, unless the G, W, or X analyte sample has a higher plate number.
# 
# Sort Replicate Filter
# The following precedence rules are applied when the analyte filter still produces more than one sample. The sort filter chooses the aliquot with the highest lexicographical sort value, to ensure that the barcode with the highest portion and/or plate number is selected when all other barcode fields are identical.
# oooooo
# Ref Link: <https://confluence.broadinstitute.org/display/GDAC/FAQ#FAQ-sampleTypesQWhatTCGAsampletypesareFirehosepipelinesexecutedupon>
#-------------------------------------------------------------------

tcga_replicateFilter = function(tsb, analyte_target=c("DNA","RNA"), decreasing=TRUE, analyte_position=20, plate=c(22,25), portion=c(18,19), filter_FFPE=FALSE, full_barcode=FALSE){
    # basically, user provide tsb and analyte_target is fine. If you
    # want to filter FFPE samples, please set filter_FFPE and full_barcode
    # all to TRUE, and tsb must have nchar of 28
    
    analyte_target = match.arg(analyte_target)
    # Strings in R are largely lexicographic
    # see ??base::Comparison
    
    # filter FFPE samples
    # provide by <http://gdac.broadinstitute.org/runs/sampleReports/latest/FPPP_FFPE_Cases.html> 
    if(full_barcode & filter_FFPE){
      ffpe = c("TCGA-44-2656-01B-06D-A271-08", "TCGA-44-2656-01B-06D-A273-01", 
               "TCGA-44-2656-01B-06D-A276-05", "TCGA-44-2656-01B-06D-A27C-26", 
               "TCGA-44-2656-01B-06R-A277-07", "TCGA-44-2662-01B-02D-A271-08", 
               "TCGA-44-2662-01B-02D-A273-01", "TCGA-44-2662-01B-02R-A277-07", 
               "TCGA-44-2665-01B-06D-A271-08", "TCGA-44-2665-01B-06D-A273-01", 
               "TCGA-44-2665-01B-06D-A276-05", "TCGA-44-2665-01B-06R-A277-07", 
               "TCGA-44-2666-01B-02D-A271-08", "TCGA-44-2666-01B-02D-A273-01", 
               "TCGA-44-2666-01B-02D-A276-05", "TCGA-44-2666-01B-02D-A27C-26", 
               "TCGA-44-2666-01B-02R-A277-07", "TCGA-44-2668-01B-02D-A271-08", 
               "TCGA-44-2668-01B-02D-A273-01", "TCGA-44-2668-01B-02D-A276-05", 
               "TCGA-44-2668-01B-02D-A27C-26", "TCGA-44-2668-01B-02R-A277-07", 
               "TCGA-44-3917-01B-02D-A271-08", "TCGA-44-3917-01B-02D-A273-01", 
               "TCGA-44-3917-01B-02D-A276-05", "TCGA-44-3917-01B-02D-A27C-26", 
               "TCGA-44-3917-01B-02R-A277-07", "TCGA-44-3918-01B-02D-A271-08", 
               "TCGA-44-3918-01B-02D-A273-01", "TCGA-44-3918-01B-02D-A276-05", 
               "TCGA-44-3918-01B-02D-A27C-26", "TCGA-44-3918-01B-02R-A277-07", 
               "TCGA-44-4112-01B-06D-A271-08", "TCGA-44-4112-01B-06D-A273-01", 
               "TCGA-44-4112-01B-06D-A276-05", "TCGA-44-4112-01B-06D-A27C-26", 
               "TCGA-44-4112-01B-06R-A277-07", "TCGA-44-5645-01B-04D-A271-08", 
               "TCGA-44-5645-01B-04D-A273-01", "TCGA-44-5645-01B-04D-A276-05", 
               "TCGA-44-5645-01B-04D-A27C-26", "TCGA-44-5645-01B-04R-A277-07", 
               "TCGA-44-6146-01B-04D-A271-08", "TCGA-44-6146-01B-04D-A273-01", 
               "TCGA-44-6146-01B-04D-A276-05", "TCGA-44-6146-01B-04D-A27C-26", 
               "TCGA-44-6146-01B-04R-A277-07", "TCGA-44-6146-01B-04R-A27D-13", 
               "TCGA-44-6147-01B-06D-A271-08", "TCGA-44-6147-01B-06D-A273-01", 
               "TCGA-44-6147-01B-06D-A276-05", "TCGA-44-6147-01B-06D-A27C-26", 
               "TCGA-44-6147-01B-06R-A277-07", "TCGA-44-6147-01B-06R-A27D-13", 
               "TCGA-44-6775-01C-02D-A271-08", "TCGA-44-6775-01C-02D-A273-01", 
               "TCGA-44-6775-01C-02D-A276-05", "TCGA-44-6775-01C-02D-A27C-26", 
               "TCGA-44-6775-01C-02R-A277-07", "TCGA-44-6775-01C-02R-A27D-13", 
               "TCGA-A6-2674-01B-04D-A270-10", "TCGA-A6-2674-01B-04R-A277-07", 
               "TCGA-A6-2677-01B-02D-A270-10", "TCGA-A6-2677-01B-02D-A274-01", 
               "TCGA-A6-2677-01B-02D-A27A-05", "TCGA-A6-2677-01B-02D-A27E-26", 
               "TCGA-A6-2677-01B-02R-A277-07", "TCGA-A6-2684-01C-08D-A270-10", 
               "TCGA-A6-2684-01C-08D-A274-01", "TCGA-A6-2684-01C-08D-A27A-05", 
               "TCGA-A6-2684-01C-08D-A27E-26", "TCGA-A6-2684-01C-08R-A277-07", 
               "TCGA-A6-3809-01B-04D-A270-10", "TCGA-A6-3809-01B-04D-A274-01", 
               "TCGA-A6-3809-01B-04D-A27A-05", "TCGA-A6-3809-01B-04D-A27E-26", 
               "TCGA-A6-3809-01B-04R-A277-07", "TCGA-A6-3810-01B-04D-A270-10", 
               "TCGA-A6-3810-01B-04D-A274-01", "TCGA-A6-3810-01B-04D-A27A-05", 
               "TCGA-A6-3810-01B-04D-A27E-26", "TCGA-A6-3810-01B-04R-A277-07", 
               "TCGA-A6-5656-01B-02D-A270-10", "TCGA-A6-5656-01B-02D-A274-01", 
               "TCGA-A6-5656-01B-02D-A27A-05", "TCGA-A6-5656-01B-02D-A27E-26", 
               "TCGA-A6-5656-01B-02R-A277-07", "TCGA-A6-5656-01B-02R-A27D-13", 
               "TCGA-A6-5659-01B-04D-A270-10", "TCGA-A6-5659-01B-04D-A274-01", 
               "TCGA-A6-5659-01B-04D-A27A-05", "TCGA-A6-5659-01B-04D-A27E-26", 
               "TCGA-A6-5659-01B-04R-A277-07", "TCGA-A6-6650-01B-02D-A270-10", 
               "TCGA-A6-6650-01B-02D-A274-01", "TCGA-A6-6650-01B-02D-A27A-05", 
               "TCGA-A6-6650-01B-02D-A27E-26", "TCGA-A6-6650-01B-02R-A277-07", 
               "TCGA-A6-6650-01B-02R-A27D-13", "TCGA-A6-6780-01B-04D-A270-10", 
               "TCGA-A6-6780-01B-04D-A274-01", "TCGA-A6-6780-01B-04D-A27A-05", 
               "TCGA-A6-6780-01B-04D-A27E-26", "TCGA-A6-6780-01B-04R-A277-07", 
               "TCGA-A6-6780-01B-04R-A27D-13", "TCGA-A6-6781-01B-06D-A270-10", 
               "TCGA-A6-6781-01B-06D-A274-01", "TCGA-A6-6781-01B-06D-A27A-05", 
               "TCGA-A6-6781-01B-06R-A277-07", "TCGA-A6-6781-01B-06R-A27D-13", 
               "TCGA-A7-A0DB-01C-02D-A272-09", "TCGA-A7-A0DB-01C-02R-A277-07", 
               "TCGA-A7-A0DB-01C-02R-A27D-13", "TCGA-A7-A13D-01B-04D-A272-09", 
               "TCGA-A7-A13D-01B-04R-A277-07", "TCGA-A7-A13D-01B-04R-A27D-13", 
               "TCGA-A7-A13E-01B-06D-A272-09", "TCGA-A7-A13E-01B-06R-A277-07", 
               "TCGA-A7-A13E-01B-06R-A27D-13", "TCGA-A7-A26E-01B-06D-A272-09", 
               "TCGA-A7-A26E-01B-06D-A275-01", "TCGA-A7-A26E-01B-06D-A27B-05", 
               "TCGA-A7-A26E-01B-06R-A277-07", "TCGA-A7-A26E-01B-06R-A27D-13", 
               "TCGA-A7-A26J-01B-02D-A272-09", "TCGA-A7-A26J-01B-02D-A275-01", 
               "TCGA-A7-A26J-01B-02D-A27B-05", "TCGA-A7-A26J-01B-02D-A27F-26", 
               "TCGA-A7-A26J-01B-02R-A277-07", "TCGA-A7-A26J-01B-02R-A27D-13", 
               "TCGA-B2-3923-01B-10D-A270-10", "TCGA-B2-3923-01B-10R-A277-07", 
               "TCGA-B2-3923-01B-10R-A27D-13", "TCGA-B2-3924-01B-03D-A270-10", 
               "TCGA-B2-3924-01B-03D-A274-01", "TCGA-B2-3924-01B-03D-A27A-05", 
               "TCGA-B2-3924-01B-03D-A27E-26", "TCGA-B2-3924-01B-03R-A277-07", 
               "TCGA-B2-3924-01B-03R-A27D-13", "TCGA-B2-5633-01B-04D-A270-10", 
               "TCGA-B2-5633-01B-04D-A274-01", "TCGA-B2-5633-01B-04D-A27A-05", 
               "TCGA-B2-5633-01B-04D-A27E-26", "TCGA-B2-5633-01B-04R-A277-07", 
               "TCGA-B2-5633-01B-04R-A27D-13", "TCGA-B2-5635-01B-04D-A270-10", 
               "TCGA-B2-5635-01B-04D-A274-01", "TCGA-B2-5635-01B-04D-A27A-05", 
               "TCGA-B2-5635-01B-04D-A27E-26", "TCGA-B2-5635-01B-04R-A277-07", 
               "TCGA-B2-5635-01B-04R-A27D-13", "TCGA-BK-A0CA-01B-02D-A272-09", 
               "TCGA-BK-A0CA-01B-02D-A275-01", "TCGA-BK-A0CA-01B-02D-A27B-05", 
               "TCGA-BK-A0CA-01B-02D-A27F-26", "TCGA-BK-A0CA-01B-02R-A277-07", 
               "TCGA-BK-A0CA-01B-02R-A27D-13", "TCGA-BK-A0CC-01B-04D-A272-09", 
               "TCGA-BK-A0CC-01B-04D-A275-01", "TCGA-BK-A0CC-01B-04D-A27B-05", 
               "TCGA-BK-A0CC-01B-04R-A277-07", "TCGA-BK-A0CC-01B-04R-A27D-13", 
               "TCGA-BK-A139-01C-08D-A272-09", "TCGA-BK-A139-01C-08D-A275-01", 
               "TCGA-BK-A139-01C-08D-A27B-05", "TCGA-BK-A139-01C-08D-A27F-26", 
               "TCGA-BK-A139-01C-08R-A277-07", "TCGA-BK-A139-01C-08R-A27D-13", 
               "TCGA-BK-A26L-01C-04D-A272-09", "TCGA-BK-A26L-01C-04D-A275-01", 
               "TCGA-BK-A26L-01C-04D-A27B-05", "TCGA-BK-A26L-01C-04D-A27F-26", 
               "TCGA-BK-A26L-01C-04R-A277-07", "TCGA-BK-A26L-01C-04R-A27D-13", 
               "TCGA-BL-A0C8-01B-04D-A271-08", "TCGA-BL-A0C8-01B-04D-A273-01", 
               "TCGA-BL-A0C8-01B-04D-A276-05", "TCGA-BL-A0C8-01B-04D-A27C-26", 
               "TCGA-BL-A0C8-01B-04R-A277-07", "TCGA-BL-A0C8-01B-04R-A27D-13", 
               "TCGA-BL-A13I-01B-04D-A271-08", "TCGA-BL-A13I-01B-04D-A276-05", 
               "TCGA-BL-A13I-01B-04R-A277-07", "TCGA-BL-A13I-01B-04R-A27D-13", 
               "TCGA-BL-A13J-01B-04D-A271-08", "TCGA-BL-A13J-01B-04D-A273-01", 
               "TCGA-BL-A13J-01B-04D-A276-05", "TCGA-BL-A13J-01B-04D-A27C-26", 
               "TCGA-BL-A13J-01B-04R-A277-07", "TCGA-BL-A13J-01B-04R-A27D-13")
      
      tsb = setdiff(tsb, tsb[which(tsb %in% ffpe)])
    }
    
    # find repeated samples
    sampleID = substr(tsb, start = 1, stop = 15)
    dp_samples = unique(sampleID[duplicated(sampleID)])
    
    if(length(dp_samples)==0){
        message("ooo Not find any duplicated barcodes, return original input..")
        tsb
    }else{
        uniq_tsb = tsb[! sampleID %in% dp_samples]
        dp_tsb = setdiff(tsb, uniq_tsb)
        
        add_tsb = c()
        
        # analyte = substr(dp_tsb, start = analyte_position, stop = analyte_position)
        # if analyte_target = "DNA"
        # analyte:  D > G,W,X
        if(analyte_target == "DNA"){
            for(x in dp_samples){
                mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
                analytes = substr(mulaliquots, 
                                  start = analyte_position,
                                  stop = analyte_position)
                if(any(analytes == "D") & !(all(analytes == "D"))){
                    aliquot = mulaliquots[which(analytes == "D")]
                    add_tsb = c(add_tsb, aliquot)
                    dp_tsb = setdiff(dp_tsb, mulaliquots)
                }
                
            }
        }else{
            for(x in dp_samples){
                mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
                analytes = substr(mulaliquots, 
                                  start = analyte_position,
                                  stop = analyte_position)
                if(any(analytes == "H") & !(all(analytes == "H"))){
                    aliquot = mulaliquots[which(analytes == "H")]
                    add_tsb = c(add_tsb, aliquot)
                    dp_tsb = setdiff(dp_tsb, mulaliquots)
                }else if(any(analytes == "R") & !(all(analytes == "R"))){
                    aliquot = mulaliquots[which(analytes == "R")]
                    add_tsb = c(add_tsb, aliquot)
                    dp_tsb = setdiff(dp_tsb, mulaliquots)
                }else if(any(analytes == "T") & !(all(analytes == "T"))){
                    aliquot = mulaliquots[which(analytes == "T")]
                    add_tsb = c(add_tsb, aliquot)
                    dp_tsb = setdiff(dp_tsb, mulaliquots)
                }
                
            }
        }
        # if analyte_target = "RNA"
        # analyte: H > R > T 
        # else{
        #     
        # }
        if(length(dp_tsb) == 0){
            message("ooo Filter barcodes successfully!")
            c(uniq_tsb, add_tsb)
        }else{
            # filter according to portion number
            sampleID_res = substr(dp_tsb, start=1, stop=15)
            dp_samples_res = unique(sampleID_res[duplicated(sampleID_res)])
            
            for(x in dp_samples_res){
                mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
                portion_codes = substr(mulaliquots,
                                       start = portion[1],
                                       stop = portion[2])
                portion_keep = sort(portion_codes, decreasing = decreasing)[1]
                if(!all(portion_codes == portion_keep)){
                    if(length(which(portion_codes == portion_keep)) == 1){
                        add_tsb = c(add_tsb, mulaliquots[which(portion_codes == portion_keep)])
                        dp_tsb = setdiff(dp_tsb, mulaliquots)
                    }else{
                        dp_tsb = setdiff(dp_tsb, mulaliquots[which(portion_codes != portion_keep)])
                    }
                    
                }
            }
            
            if(length(dp_tsb)==0){
                message("ooo Filter barcodes successfully!")
                c(uniq_tsb, add_tsb)
            }else{
                # filter according to plate number
                sampleID_res = substr(dp_tsb, start=1, stop=15)
                dp_samples_res = unique(sampleID_res[duplicated(sampleID_res)])
                for(x in dp_samples_res){
                    mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
                    plate_codes = substr(mulaliquots,
                                         start = plate[1],
                                         stop = plate[2])
                    plate_keep = sort(plate_codes, decreasing = decreasing)[1]
                    add_tsb = c(add_tsb, mulaliquots[which(plate_codes == plate_keep)])
                    dp_tsb = setdiff(dp_tsb, mulaliquots)
                }
                
                if(length(dp_tsb)==0){
                    message("ooo Filter barcodes successfully!")
                    c(uniq_tsb, add_tsb)
                }else{
                    message("ooo barcodes ", dp_tsb, " failed in filter process, other barcodes will be returned.")
                    c(uniq_tsb, add_tsb)
                }
            }
        }
    }
}


#-------------------------
# Test and Check result 
#------------------------

tsb = c("TCGA-55-7995-01A-11D-2183-01", "TCGA-55-8616-01A-11D-2389-01", 
        "TCGA-NJ-A4YG-01A-22D-A25M-01", "TCGA-44-2668-01A-01D-A273-01", 
        "TCGA-49-6761-01A-31D-1943-01", "TCGA-75-6207-01A-11D-1752-01", 
        "TCGA-73-4659-01A-01D-1204-01", "TCGA-78-7536-01A-11D-2062-01", 
        "TCGA-J2-8192-01A-11D-2237-01", "TCGA-49-AAR9-01A-21D-A40Z-01", 
        "TCGA-L4-A4E5-01A-11D-A24O-01", "TCGA-95-A4VK-01A-11D-A25M-01", 
        "TCGA-49-4486-01A-01D-1204-01", "TCGA-67-3772-01A-01D-0944-01", 
        "TCGA-78-7159-01A-11D-2035-01", "TCGA-55-A491-01A-11D-A24C-01", 
        "TCGA-05-4403-01A-01D-1204-01", "TCGA-97-8171-01A-11D-2283-01", 
        "TCGA-L4-A4E6-01A-11D-A24C-01", "TCGA-49-6743-01A-11D-1854-01", 
        "TCGA-62-8397-01A-11D-2322-01", "TCGA-MP-A4TF-01A-11D-A25M-01", 
        "TCGA-55-8207-01A-11D-2237-01", "TCGA-53-7813-01A-11D-2166-01", 
        "TCGA-44-7671-01A-11D-2062-01", "TCGA-78-7145-01A-11D-2035-01", 
        "TCGA-78-8655-01A-11D-2389-01", "TCGA-86-7955-01A-11D-2183-01", 
        "TCGA-4B-A93V-01A-11D-A396-01", "TCGA-MP-A4T6-01A-32D-A25M-01", 
        "TCGA-55-6986-01A-11D-1943-01", "TCGA-55-6981-01A-11D-1943-01", 
        "TCGA-97-8179-01A-11D-2283-01", "TCGA-44-5645-01A-01D-1624-01", 
        "TCGA-78-7537-01A-11D-2062-01", "TCGA-69-7760-01A-11D-2166-01", 
        "TCGA-55-6975-01A-11D-1943-01", "TCGA-49-6744-01A-11D-1854-01", 
        "TCGA-44-6776-01A-11D-1854-01", "TCGA-73-4677-01A-01D-1204-01", 
        "TCGA-75-7031-01A-11D-1943-01", "TCGA-55-7815-01A-11D-2166-01", 
        "TCGA-J2-A4AE-01A-21D-A24C-01", "TCGA-69-8253-01A-11D-2283-01", 
        "TCGA-38-6178-01A-11D-1752-01", "TCGA-55-8514-01A-11D-2389-01", 
        "TCGA-55-8508-01A-11D-2389-01", "TCGA-91-6828-01A-11D-1854-01", 
        "TCGA-78-7149-01A-11D-2035-01", "TCGA-64-5779-01A-01D-1624-01", 
        "TCGA-62-8399-01A-21D-2322-01", "TCGA-44-6779-01A-11D-1854-01", 
        "TCGA-49-6745-01A-11D-1854-01", "TCGA-44-5643-01A-01D-1624-01", 
        "TCGA-62-A46U-01A-11D-A24C-01", "TCGA-55-8506-01A-11D-2389-01", 
        "TCGA-44-A4SU-01A-11D-A24O-01", "TCGA-55-8301-01A-11D-2283-01", 
        "TCGA-05-5425-01A-02D-1624-01", "TCGA-55-8619-01A-11D-2389-01", 
        "TCGA-86-A4P7-01A-11D-A24O-01", "TCGA-55-8208-01A-11D-2237-01", 
        "TCGA-05-4398-01A-01D-1204-01", "TCGA-55-7728-01A-11D-2183-01", 
        "TCGA-49-AAR0-01A-21D-A396-01", "TCGA-44-3398-01A-01D-1549-01", 
        "TCGA-78-7540-01A-11D-2062-01", "TCGA-64-5815-01A-01D-1624-01", 
        "TCGA-62-A46Y-01A-11D-A24C-01", "TCGA-53-7624-01A-11D-2062-01", 
        "TCGA-73-4676-01A-01D-1752-01", "TCGA-62-8395-01A-11D-2322-01", 
        "TCGA-55-6985-01A-11D-1943-01", "TCGA-49-AARO-01A-12D-A40Z-01", 
        "TCGA-95-7567-01A-11D-2062-01", "TCGA-75-5126-01A-01D-1752-01", 
        "TCGA-86-7954-01A-11D-2183-01", "TCGA-86-7701-01A-11D-2166-01", 
        "TCGA-NJ-A4YP-01A-11D-A25M-01", "TCGA-95-7562-01A-11D-2237-01", 
        "TCGA-44-2655-01A-01D-0944-01", "TCGA-05-4390-01A-02D-1752-01", 
        "TCGA-NJ-A4YI-01A-11D-A25M-01", "TCGA-50-8457-01A-11D-2322-01", 
        "TCGA-MP-A4TI-01A-21D-A24O-01", "TCGA-97-A4M7-01A-11D-A24O-01", 
        "TCGA-78-7155-01A-11D-2035-01", "TCGA-71-6725-01A-11D-1854-01", 
        "TCGA-55-7726-01A-11D-2166-01", "TCGA-50-6595-01A-12D-1854-01", 
        "TCGA-44-2657-01A-01D-1877-01", "TCGA-78-7158-01A-11D-2035-01", 
        "TCGA-86-8672-01A-21D-2389-01", "TCGA-55-8507-01A-11D-2389-01", 
        "TCGA-75-6206-01A-11D-1752-01", "TCGA-MP-A4SY-01A-21D-A24O-01", 
        "TCGA-NJ-A7XG-01A-12D-A396-01", "TCGA-05-4422-01A-01D-1204-01", 
        "TCGA-97-8547-01A-11D-2389-01", "TCGA-38-4631-01A-01D-1752-01", 
        "TCGA-75-6203-01A-11D-1752-01", "TCGA-97-7937-01A-11D-2166-01", 
        "TCGA-05-4396-01A-21D-1854-01", "TCGA-97-A4M5-01A-11D-A24O-01", 
        "TCGA-62-8398-01A-11D-2322-01", "TCGA-55-8094-01A-11D-2237-01", 
        "TCGA-86-8281-01A-11D-2283-01", "TCGA-MP-A4TD-01A-32D-A25M-01", 
        "TCGA-44-2655-01A-01D-1549-01", "TCGA-44-2668-01B-02D-A273-01", 
        "TCGA-55-6712-01A-11D-1854-01", "TCGA-55-8092-01A-11D-2237-01", 
        "TCGA-49-4494-01A-01D-1204-01", "TCGA-91-6830-01A-11D-1943-01", 
        "TCGA-05-4425-01A-01D-1752-01", "TCGA-55-8511-01A-11D-2389-01", 
        "TCGA-49-4501-01A-01D-1204-01", "TCGA-86-8668-01A-11D-2389-01", 
        "TCGA-NJ-A55O-01A-11D-A25M-01", "TCGA-86-A4D0-01A-11D-A24C-01", 
        "TCGA-S2-AA1A-01A-12D-A396-01", "TCGA-64-1676-01A-01D-0944-01", 
        "TCGA-49-4490-01A-21D-1854-01", "TCGA-69-7761-01A-11D-2166-01", 
        "TCGA-50-5051-01A-21D-1854-01", "TCGA-86-8669-01A-11D-2389-01", 
        "TCGA-62-A472-01A-11D-A24C-01", "TCGA-44-6775-01C-02D-A273-01", 
        "TCGA-44-2665-01A-01D-0944-01", "TCGA-44-6147-01B-06D-A273-01", 
        "TCGA-44-3919-01A-02D-1877-01", "TCGA-75-5125-01A-01D-1752-01", 
        "TCGA-80-5611-01A-01D-1624-01", "TCGA-50-5939-01A-11D-1624-01", 
        "TCGA-69-7765-01A-11D-2166-01", "TCGA-86-7713-01A-11D-2062-01", 
        "TCGA-05-5715-01A-01D-1624-01", "TCGA-MP-A4TC-01A-11D-A24O-01", 
        "TCGA-49-4505-01A-01D-1204-01", "TCGA-49-6742-01A-11D-1854-01", 
        "TCGA-99-8028-01A-11D-2237-01", "TCGA-38-4632-01A-01D-1752-01", 
        "TCGA-67-6217-01A-11D-1752-01", "TCGA-55-6970-01A-11D-1943-01", 
        "TCGA-64-1680-01A-02D-0944-01", "TCGA-73-4670-01A-01D-1204-01", 
        "TCGA-50-5936-01A-11D-1624-01", "TCGA-55-7725-01A-11D-2166-01", 
        "TCGA-49-4514-01A-21D-1854-01", "TCGA-50-5068-01A-01D-1624-01", 
        "TCGA-53-A4EZ-01A-12D-A24O-01", "TCGA-05-4426-01A-01D-1204-01", 
        "TCGA-44-2656-01A-02D-1549-01", "TCGA-55-A493-01A-11D-A24C-01", 
        "TCGA-05-5428-01A-01D-1624-01", "TCGA-44-2659-01A-01D-0944-01", 
        "TCGA-44-2662-01A-01D-0944-01", "TCGA-69-8255-01A-11D-2283-01", 
        "TCGA-38-4630-01A-01D-1204-01", "TCGA-44-2659-01A-01D-1549-01", 
        "TCGA-69-7764-01A-11D-2166-01", "TCGA-75-7027-01A-11D-1943-01", 
        "TCGA-55-6642-01A-11D-1854-01", "TCGA-75-7030-01A-11D-1943-01", 
        "TCGA-97-8177-01A-11D-2283-01", "TCGA-MP-A4T8-01A-11D-A24O-01", 
        "TCGA-50-6673-01A-11D-1943-01", "TCGA-44-7669-01A-21D-2062-01", 
        "TCGA-05-4395-01A-01D-1204-01", "TCGA-44-4112-01B-06D-A273-01", 
        "TCGA-99-8033-01A-11D-2237-01", "TCGA-44-6146-01B-04D-A273-01", 
        "TCGA-38-4627-01A-01D-1549-01", "TCGA-55-8620-01A-11D-2389-01", 
        "TCGA-05-4397-01A-01D-1204-01", "TCGA-50-5941-01A-11D-1752-01", 
        "TCGA-73-7499-01A-11D-2183-01", "TCGA-NJ-A55A-01A-11D-A25M-01", 
        "TCGA-67-3774-01A-01D-0944-01", "TCGA-95-7948-01A-11D-2183-01", 
        "TCGA-95-8494-01A-11D-2322-01", "TCGA-62-A46O-01A-11D-A24C-01", 
        "TCGA-44-6147-01A-11D-A273-01", "TCGA-44-2662-01A-01D-A273-01", 
        "TCGA-50-5944-01A-11D-1752-01", "TCGA-86-8671-01A-11D-2389-01", 
        "TCGA-75-6214-01A-41D-1943-01", "TCGA-55-8085-01A-11D-2237-01", 
        "TCGA-44-7660-01A-11D-2062-01", "TCGA-49-4488-01A-01D-1752-01", 
        "TCGA-55-6969-01A-11D-1943-01", "TCGA-97-8174-01A-11D-2283-01", 
        "TCGA-75-5122-01A-01D-1752-01", "TCGA-44-7667-01A-31D-2062-01", 
        "TCGA-91-6831-01A-11D-1854-01", "TCGA-64-1677-01A-01D-0944-01", 
        "TCGA-44-3918-01A-01D-1877-01", "TCGA-50-5933-01A-11D-1752-01", 
        "TCGA-55-8505-01A-11D-2389-01", "TCGA-50-5946-01A-11D-1752-01", 
        "TCGA-75-6212-01A-11D-1752-01", "TCGA-78-8660-01A-11D-2389-01", 
        "TCGA-97-7553-01A-21D-2035-01", "TCGA-05-4420-01A-01D-1204-01", 
        "TCGA-44-A47F-01A-11D-A24C-01", "TCGA-50-6591-01A-11D-1752-01", 
        "TCGA-38-A44F-01A-11D-A24C-01", "TCGA-44-5645-01B-04D-A273-01", 
        "TCGA-44-3396-01A-01D-1549-01", "TCGA-44-6775-01A-11D-1854-01", 
        "TCGA-44-2662-01B-02D-A273-01", "TCGA-55-6972-01A-11D-1943-01", 
        "TCGA-50-5045-01A-01D-1624-01", "TCGA-55-6971-01A-11D-1943-01", 
        "TCGA-55-A57B-01A-12D-A396-01", "TCGA-55-8087-01A-11D-2237-01", 
        "TCGA-78-7163-01A-12D-2062-01", "TCGA-62-A46P-01A-11D-A24C-01", 
        "TCGA-MP-A4T2-01A-11D-A24O-01", "TCGA-44-A479-01A-31D-A24C-01", 
        "TCGA-MP-A4SW-01A-21D-A24O-01", "TCGA-05-4415-01A-22D-1854-01", 
        "TCGA-53-7626-01A-12D-2062-01", "TCGA-55-7727-01A-11D-2166-01", 
        "TCGA-55-6984-01A-11D-1943-01", "TCGA-MP-A4T4-01A-11D-A25M-01", 
        "TCGA-44-6775-01A-11D-A273-01", "TCGA-55-1594-01A-01D-0944-01", 
        "TCGA-78-7147-01A-11D-2035-01", "TCGA-67-3770-01A-01D-0944-01", 
        "TCGA-NJ-A4YF-01A-12D-A25M-01", "TCGA-73-4666-01A-01D-1204-01", 
        "TCGA-97-A4M6-01A-11D-A24O-01", "TCGA-86-8076-01A-31D-2237-01", 
        "TCGA-49-AARR-01A-11D-A40Z-01", "TCGA-44-2666-01B-02D-A273-01", 
        "TCGA-55-8621-01A-11D-2389-01", "TCGA-78-7539-01A-11D-2062-01", 
        "TCGA-50-8459-01A-11D-2322-01", "TCGA-78-7160-01A-11D-2035-01", 
        "TCGA-35-5375-01A-01D-1624-01", "TCGA-38-4625-01A-01D-1549-01", 
        "TCGA-78-7148-01A-11D-2035-01", "TCGA-05-4417-01A-22D-1854-01", 
        "TCGA-50-8460-01A-11D-2322-01", "TCGA-J2-A4AG-01A-11D-A24C-01", 
        "TCGA-55-8204-01A-11D-2237-01", "TCGA-67-3773-01A-01D-0944-01", 
        "TCGA-44-2666-01A-01D-A273-01", "TCGA-NJ-A4YQ-01A-11D-A25M-01", 
        "TCGA-86-8674-01A-21D-2389-01", "TCGA-91-6849-01A-11D-1943-01", 
        "TCGA-75-6211-01A-11D-1752-01", "TCGA-38-4626-01A-01D-1204-01", 
        "TCGA-44-7661-01A-11D-2062-01", "TCGA-55-7284-01B-11D-2237-01", 
        "TCGA-38-4626-01A-01D-1549-01", "TCGA-97-7554-01A-11D-2035-01", 
        "TCGA-49-AAR3-01A-11D-A40Z-01", "TCGA-69-A59K-01A-11D-A25M-01", 
        "TCGA-50-6592-01A-11D-1752-01", "TCGA-55-A4DG-01A-11D-A24C-01", 
        "TCGA-99-AA5R-01A-11D-A396-01", "TCGA-86-8585-01A-11D-2389-01", 
        "TCGA-67-4679-01B-01D-1752-01", "TCGA-55-8203-01A-11D-2237-01", 
        "TCGA-95-7039-01A-11D-1943-01", "TCGA-MP-A4TE-01A-22D-A25M-01", 
        "TCGA-97-A4M2-01A-12D-A24O-01", "TCGA-O1-A52J-01A-11D-A25M-01", 
        "TCGA-78-7167-01A-11D-2062-01", "TCGA-55-6980-01A-11D-1943-01", 
        "TCGA-86-8358-01A-11D-2322-01", "TCGA-67-6216-01A-11D-1752-01", 
        "TCGA-44-3918-01A-01D-A273-01", "TCGA-69-8453-01A-12D-2322-01", 
        "TCGA-05-4244-01A-01D-1877-01", "TCGA-93-A4JO-01A-21D-A24O-01", 
        "TCGA-50-5931-01A-11D-1752-01", "TCGA-55-8299-01A-11D-2283-01", 
        "TCGA-MP-A4TA-01A-21D-A24O-01", "TCGA-67-6215-01A-11D-1752-01", 
        "TCGA-97-7941-01A-11D-2183-01", "TCGA-50-5055-01A-01D-1624-01", 
        "TCGA-93-8067-01A-11D-2283-01", "TCGA-55-7576-01A-11D-2062-01", 
        "TCGA-44-2657-01A-01D-1549-01", "TCGA-NJ-A55R-01A-11D-A25M-01", 
        "TCGA-78-7150-01A-21D-2035-01", "TCGA-69-7763-01A-11D-2166-01", 
        "TCGA-44-4112-01A-01D-A273-01", "TCGA-55-8090-01A-11D-2237-01", 
        "TCGA-69-7979-01A-11D-2183-01", "TCGA-55-5899-01A-11D-1624-01", 
        "TCGA-97-8552-01A-11D-2389-01", "TCGA-05-4432-01A-01D-1204-01", 
        "TCGA-86-8074-01A-11D-2237-01", "TCGA-05-4382-01A-01D-1204-01", 
        "TCGA-55-8089-01A-11D-2237-01", "TCGA-78-7153-01A-11D-2035-01", 
        "TCGA-44-2662-01A-01D-1549-01", "TCGA-55-7910-01A-11D-2166-01", 
        "TCGA-67-3771-01A-01D-0944-01", "TCGA-55-7574-01A-11D-2035-01", 
        "TCGA-49-AARQ-01A-11D-A40Z-01", "TCGA-44-6777-01A-11D-1854-01", 
        "TCGA-97-7552-01A-11D-2035-01", "TCGA-44-2656-01A-02D-0944-01", 
        "TCGA-86-8055-01A-11D-2237-01", "TCGA-86-7714-01A-12D-2166-01", 
        "TCGA-44-2668-01A-01D-0944-01", "TCGA-L9-A743-01A-43D-A396-01", 
        "TCGA-62-A470-01A-11D-A24C-01", "TCGA-97-7938-01A-11D-2166-01", 
        "TCGA-05-5429-01A-01D-1624-01", "TCGA-44-6145-01A-11D-1752-01", 
        "TCGA-05-5423-01A-01D-1624-01", "TCGA-49-4510-01A-01D-1204-01", 
        "TCGA-71-8520-01A-11D-2389-01", "TCGA-62-A46S-01A-11D-A24C-01", 
        "TCGA-55-6987-01A-11D-1943-01", "TCGA-55-6983-01A-11D-1943-01", 
        "TCGA-55-6968-01A-11D-1943-01", "TCGA-44-2661-01A-01D-1877-01", 
        "TCGA-97-7546-01A-11D-2035-01", "TCGA-64-5778-01A-01D-1624-01", 
        "TCGA-44-2668-01A-01D-1549-01", "TCGA-50-5066-01A-01D-1624-01", 
        "TCGA-97-A4M1-01A-11D-A24O-01", "TCGA-86-8073-01A-11D-2237-01", 
        "TCGA-86-8280-01A-11D-2283-01", "TCGA-91-6847-01A-11D-1943-01", 
        "TCGA-49-4487-01A-21D-1854-01", "TCGA-L9-A444-01A-21D-A24C-01", 
        "TCGA-80-5608-01A-31D-1943-01", "TCGA-75-6205-01A-11D-1752-01", 
        "TCGA-05-4427-01A-21D-1854-01", "TCGA-49-AARN-01A-21D-A40Z-01", 
        "TCGA-49-6767-01A-11D-1854-01", "TCGA-78-7161-01A-11D-2035-01", 
        "TCGA-MP-A4TK-01A-11D-A24O-01", "TCGA-05-4389-01A-01D-1204-01", 
        "TCGA-64-5775-01A-01D-1624-01", "TCGA-MP-A5C7-01A-11D-A25M-01", 
        "TCGA-44-5645-01A-01D-A273-01", "TCGA-05-4249-01A-01D-1877-01", 
        "TCGA-49-AARE-01A-11D-A40Z-01", "TCGA-75-5147-01A-01D-1624-01", 
        "TCGA-78-7156-01A-11D-2035-01", "TCGA-97-8172-01A-11D-2283-01", 
        "TCGA-95-7944-01A-11D-2183-01", "TCGA-44-3917-01A-01D-A273-01", 
        "TCGA-78-7535-01A-11D-2062-01", "TCGA-38-4625-01A-01D-1204-01", 
        "TCGA-91-8499-01A-11D-2389-01", "TCGA-05-4405-01A-21D-1854-01", 
        "TCGA-78-8640-01A-11D-2389-01", "TCGA-93-A4JP-01A-11D-A24O-01", 
        "TCGA-55-7573-01A-11D-2035-01", "TCGA-91-A4BC-01A-11D-A24C-01", 
        "TCGA-44-3396-01A-01D-1204-01", "TCGA-L9-A443-01A-12D-A24C-01", 
        "TCGA-44-2665-01A-01D-1549-01", "TCGA-50-5049-01A-01D-1624-01", 
        "TCGA-44-6146-01A-11D-1752-01", "TCGA-95-7043-01A-11D-1943-01", 
        "TCGA-64-5781-01A-01D-1624-01", "TCGA-MN-A4N1-01A-11D-A24O-01", 
        "TCGA-69-7973-01A-11D-2183-01", "TCGA-55-8206-01A-11D-2237-01", 
        "TCGA-86-8673-01A-11D-2389-01", "TCGA-MN-A4N4-01A-12D-A24O-01", 
        "TCGA-55-7914-01A-11D-2166-01", "TCGA-91-6840-01A-11D-1943-01", 
        "TCGA-05-4433-01A-22D-1854-01", "TCGA-55-8512-01A-11D-2389-01", 
        "TCGA-44-2656-01B-06D-A273-01", "TCGA-44-2665-01A-01D-A273-01", 
        "TCGA-35-4123-01A-01D-1877-01", "TCGA-91-8496-01A-11D-2389-01", 
        "TCGA-78-7152-01A-11D-2035-01", "TCGA-50-6594-01A-11D-1752-01", 
        "TCGA-49-4506-01A-01D-1204-01", "TCGA-83-5908-01A-21D-2283-01", 
        "TCGA-97-A4LX-01A-11D-A24O-01", "TCGA-MP-A4SV-01A-11D-A24O-01", 
        "TCGA-55-A48Y-01A-11D-A24C-01", "TCGA-44-A47G-01A-21D-A24C-01", 
        "TCGA-55-7724-01A-11D-2166-01", "TCGA-78-7220-01A-11D-2035-01", 
        "TCGA-86-A456-01A-11D-A24C-01", "TCGA-91-8497-01A-11D-2389-01", 
        "TCGA-55-A494-01A-11D-A24O-01", "TCGA-05-4384-01A-01D-1752-01", 
        "TCGA-44-3398-01A-01D-1877-01", "TCGA-62-8402-01A-11D-2322-01", 
        "TCGA-38-4627-01A-01D-1204-01", "TCGA-69-7980-01A-11D-2183-01", 
        "TCGA-55-7570-01A-11D-2035-01", "TCGA-55-A48X-01A-11D-A24C-01", 
        "TCGA-05-4418-01A-01D-1204-01", "TCGA-55-7283-01A-11D-2035-01", 
        "TCGA-91-6829-01A-21D-1854-01", "TCGA-86-8075-01A-11D-2237-01", 
        "TCGA-55-6978-01A-11D-1943-01", "TCGA-64-1678-01A-01D-0944-01", 
        "TCGA-J2-8194-01A-11D-2237-01", "TCGA-05-4402-01A-01D-1204-01", 
        "TCGA-50-5932-01A-11D-1752-01", "TCGA-44-2661-01A-01D-1549-01", 
        "TCGA-49-4507-01A-01D-1204-01", "TCGA-97-7547-01A-11D-2035-01", 
        "TCGA-86-7953-01A-11D-2183-01", "TCGA-55-8096-01A-11D-2237-01", 
        "TCGA-95-8039-01A-11D-2237-01", "TCGA-55-8510-01A-11D-2389-01", 
        "TCGA-93-A4JN-01A-11D-A24O-01", "TCGA-86-A4JF-01A-11D-A24O-01", 
        "TCGA-05-5420-01A-01D-1624-01", "TCGA-55-8302-01A-11D-2322-01", 
        "TCGA-99-8032-01A-11D-2237-01", "TCGA-78-7166-01A-12D-2062-01", 
        "TCGA-55-7907-01A-11D-2166-01", "TCGA-05-4250-01A-01D-1877-01", 
        "TCGA-50-5930-01A-11D-1752-01", "TCGA-55-8614-01A-11D-2389-01", 
        "TCGA-38-4629-01A-02D-1204-01", "TCGA-MN-A4N5-01A-11D-A24O-01", 
        "TCGA-62-A46V-01A-11D-A24C-01", "TCGA-55-6543-01A-11D-1752-01", 
        "TCGA-69-7974-01A-11D-2183-01", "TCGA-44-A4SS-01A-11D-A24O-01", 
        "TCGA-78-7633-01A-11D-2062-01", "TCGA-44-A47A-01A-21D-A24C-01", 
        "TCGA-MP-A4T7-01A-11D-A24O-01", "TCGA-55-A492-01A-11D-A24C-01", 
        "TCGA-86-8056-01A-11D-2237-01", "TCGA-86-8278-01A-11D-2283-01", 
        "TCGA-35-3615-01A-01D-0944-01", "TCGA-44-6774-01A-21D-1854-01", 
        "TCGA-35-4122-01A-01D-1877-01", "TCGA-78-8648-01A-11D-2389-01", 
        "TCGA-L9-A50W-01A-12D-A396-01", "TCGA-86-A4P8-01A-11D-A24O-01", 
        "TCGA-91-6836-01A-21D-1854-01", "TCGA-05-4430-01A-02D-1204-01", 
        "TCGA-93-7347-01A-11D-2183-01", "TCGA-64-5774-01A-01D-1624-01", 
        "TCGA-55-1595-01A-01D-0944-01", "TCGA-44-8119-01A-11D-2237-01", 
        "TCGA-73-7498-01A-12D-2183-01", "TCGA-50-5942-01A-21D-1752-01", 
        "TCGA-50-6597-01A-11D-1854-01", "TCGA-44-7670-01A-11D-2062-01", 
        "TCGA-MP-A4T9-01A-11D-A24O-01", "TCGA-78-7143-01A-11D-2035-01", 
        "TCGA-44-7659-01A-11D-2062-01", "TCGA-75-7025-01A-12D-1943-01", 
        "TCGA-44-6147-01A-11D-1752-01", "TCGA-62-8394-01A-11D-2322-01", 
        "TCGA-80-5607-01A-31D-1943-01", "TCGA-44-8120-01A-11D-2237-01", 
        "TCGA-69-7978-01A-11D-2183-01", "TCGA-91-7771-01A-11D-2166-01", 
        "TCGA-44-2665-01B-06D-A273-01", "TCGA-86-8279-01A-11D-2283-01", 
        "TCGA-97-A4M3-01A-11D-A24O-01", "TCGA-73-4658-01A-01D-1752-01", 
        "TCGA-55-1592-01A-01D-0944-01", "TCGA-55-8205-01A-11D-2237-01", 
        "TCGA-86-6851-01A-11D-1943-01", "TCGA-55-8097-01A-11D-2237-01", 
        "TCGA-86-8359-01A-11D-2322-01", "TCGA-64-1679-01A-21D-2062-01", 
        "TCGA-86-7711-01A-11D-2062-01", "TCGA-55-7281-01A-11D-2035-01", 
        "TCGA-44-2666-01A-01D-0944-01", "TCGA-44-2656-01A-02D-A273-01", 
        "TCGA-55-8615-01A-11D-2389-01", "TCGA-J2-A4AD-01A-11D-A24C-01", 
        "TCGA-75-5146-01A-01D-1624-01", "TCGA-38-7271-01A-11D-2035-01", 
        "TCGA-55-7227-01A-11D-2035-01", "TCGA-93-A4JQ-01A-11D-A24O-01", 
        "TCGA-44-7672-01A-11D-2062-01", "TCGA-55-6979-01A-11D-1943-01", 
        "TCGA-86-6562-01A-11D-1752-01", "TCGA-MP-A4TJ-01A-51D-A25M-01", 
        "TCGA-99-8025-01A-11D-2237-01", "TCGA-38-4628-01A-01D-1204-01", 
        "TCGA-05-4424-01A-22D-1854-01", "TCGA-97-8176-01A-11D-2389-01", 
        "TCGA-97-8175-01A-11D-2283-01", "TCGA-55-7911-01A-11D-2166-01", 
        "TCGA-05-4434-01A-01D-1204-01", "TCGA-44-6146-01A-11D-A273-01", 
        "TCGA-73-4668-01A-01D-1204-01", "TCGA-86-8054-01A-11D-2237-01", 
        "TCGA-50-5044-01A-21D-1854-01", "TCGA-44-6778-01A-11D-1854-01", 
        "TCGA-62-A46R-01A-11D-A24C-01", "TCGA-44-3917-01B-02D-A273-01", 
        "TCGA-78-7542-01A-21D-2062-01", "TCGA-64-1681-01A-11D-2062-01", 
        "TCGA-91-6848-01A-11D-1943-01", "TCGA-55-8513-01A-11D-2389-01", 
        "TCGA-L9-A7SV-01A-11D-A396-01", "TCGA-49-AAQV-01A-11D-A396-01", 
        "TCGA-49-AAR2-01A-11D-A396-01", "TCGA-44-4112-01A-01D-1877-01", 
        "TCGA-05-4410-01A-21D-1854-01", "TCGA-MP-A4TH-01A-31D-A25M-01", 
        "TCGA-55-A4DF-01A-11D-A24C-01", "TCGA-78-7162-01A-21D-2062-01", 
        "TCGA-44-5644-01A-21D-2035-01", "TCGA-73-A9RS-01A-11D-A40Z-01", 
        "TCGA-55-A490-01A-11D-A24C-01", "TCGA-55-8091-01A-11D-2237-01", 
        "TCGA-49-AAR4-01A-12D-A40Z-01", "TCGA-78-7146-01A-11D-2035-01", 
        "TCGA-95-A4VP-01A-21D-A25M-01", "TCGA-93-7348-01A-21D-2035-01", 
        "TCGA-69-8254-01A-11D-2283-01", "TCGA-44-7662-01A-11D-2062-01", 
        "TCGA-78-8662-01A-11D-2389-01", "TCGA-44-3918-01B-02D-A273-01", 
        "TCGA-91-A4BD-01A-11D-A24C-01", "TCGA-95-7947-01A-11D-2183-01", 
        "TCGA-55-A48Z-01A-12D-A24O-01", "TCGA-49-4512-01A-21D-1854-01", 
        "TCGA-73-4662-01A-01D-1204-01", "TCGA-50-5072-01A-21D-1854-01", 
        "TCGA-95-A4VN-01A-11D-A25M-01", "TCGA-55-1596-01A-01D-0944-01", 
        "TCGA-44-A47B-01A-11D-A24C-01", "TCGA-97-A4M0-01A-11D-A24O-01", 
        "TCGA-44-6148-01A-11D-1752-01", "TCGA-55-7903-01A-11D-2166-01", 
        "TCGA-50-6593-01A-11D-1752-01", "TCGA-L9-A8F4-01A-11D-A396-01", 
        "TCGA-55-7994-01A-11D-2183-01", "TCGA-L9-A5IP-01A-21D-A396-01", 
        "TCGA-55-6982-01A-11D-1943-01", "TCGA-44-8117-01A-11D-2237-01", 
        "TCGA-50-5935-01A-11D-1752-01", "TCGA-78-7154-01A-11D-2035-01", 
        "TCGA-99-7458-01A-11D-2035-01", "TCGA-62-A471-01A-12D-A24C-01", 
        "TCGA-50-6590-01A-12D-1854-01", "TCGA-91-6835-01A-11D-1854-01", 
        "TCGA-44-6144-01A-11D-1752-01", "TCGA-73-4675-01A-01D-1204-01", 
        "TCGA-50-7109-01A-11D-2035-01", "TCGA-55-7913-01B-11D-2237-01"
)


# Sort replicate filter
require(tidyverse)
tsb %>% substr(1,15) %>% unique() -> test_tsb
    
uniq_tsb = tcga_replicateFilter(tsb = tsb)

all(test_tsb %in% substr(uniq_tsb, 1, 15))
all.equal(sort(test_tsb), sort(substr(uniq_tsb, 1, 15)))

# DNA analyte replicate filter
tcga_replicateFilter(tsb = c("TCGA-55-7913-01B-11D-2237-01", "TCGA-55-7913-01B-11X-2237-01", "TCGA-55-7913-01B-11D-2237-01"))

# RNA analyte replicate filter
tcga_replicateFilter(tsb = c("TCGA-55-7913-01B-11H-2237-01", "TCGA-55-7913-01B-11R-2237-01", "TCGA-55-7913-01B-11T-2237-01"), analyte_target = "RNA")
tcga_replicateFilter(tsb = c("TCGA-55-7913-01B-11R-2237-01", "TCGA-55-7913-01B-11R-2237-01", "TCGA-55-7913-01B-11T-2237-01"), analyte_target = "RNA")
tcga_replicateFilter(tsb = c("TCGA-55-7913-01B-11T-2237-01", "TCGA-55-7913-01B-11T-2237-01", "TCGA-55-7913-01B-11D-2237-01"), analyte_target = "RNA")




####### test result
#> all(test_tsb %in% substr(uniq_tsb, 1, 15))
#[1] TRUE
#> all.equal(sort(test_tsb), sort(substr(uniq_tsb, 1, 15)))
#[1] TRUE
#> # DNA analyte replicate filter
#> tcga_replicateFilter(tsb = c("TCGA-55-7913-01B-11D-2237-01", "TCGA-55-7913-01B-11X-2237-01", "TCGA-55-7913-01B-11D-2237-01"))
#ooo Filter barcodes successfully!
#[1] "TCGA-55-7913-01B-11D-2237-01"
#> # RNA analyte replicate filter
#> tcga_replicateFilter(tsb = c("TCGA-55-7913-01B-11H-2237-01", "TCGA-55-7913-01B-11R-2237-01", "TCGA-55-7913-01B-11T-2237-01"), analyte_target = "RNA")
#ooo Filter barcodes successfully!
#[1] "TCGA-55-7913-01B-11H-2237-01"
#> tcga_replicateFilter(tsb = c("TCGA-55-7913-01B-11R-2237-01", "TCGA-55-7913-01B-11R-2237-01", "TCGA-55-7913-01B-11T-2237-01"), analyte_target = "RNA")
#ooo Filter barcodes successfully!
#[1] "TCGA-55-7913-01B-11R-2237-01"
#> tcga_replicateFilter(tsb = c("TCGA-55-7913-01B-11T-2237-01", "TCGA-55-7913-01B-11T-2237-01", "TCGA-55-7913-01B-11D-2237-01"), analyte_target = "RNA")
#ooo Filter barcodes successfully!
#[1] "TCGA-55-7913-01B-11T-2237-01"

## filter FFPE samples, provied in <http://gdac.broadinstitute.org/runs/sampleReports/latest/FPPP_FFPE_Cases.html>
## please note 
    # basically, user provide tsb and analyte_target is fine. If you
    # want to filter FFPE samples, please set filter_FFPE and full_barcode
    # all to TRUE, and tsb must have nchar of 28
# > tcga_replicateFilter(tsb = c("TCGA-55-7913-01B-11H-2237-01", "TCGA-55-7913-01B-11R-2237-01", "TCGA-55-7913-01B-11T-2237-01", "TCGA-44-2656-01B-06D-A271-08", "TCGA-44-2656-01B-06D-A273-01"), analyte_target = "RNA")
# ooo Filter barcodes successfully!
#   [1] "TCGA-55-7913-01B-11H-2237-01" "TCGA-44-2656-01B-06D-A273-01"
# > tcga_replicateFilter(tsb = c("TCGA-55-7913-01B-11H-2237-01", "TCGA-55-7913-01B-11R-2237-01", "TCGA-55-7913-01B-11T-2237-01", "TCGA-44-2656-01B-06D-A271-08", "TCGA-44-2656-01B-06D-A273-01"), analyte_target = "RNA", filter_FFPE = TRUE, full_barcode = TRUE)
# ooo Filter barcodes successfully!
#   [1] "TCGA-55-7913-01B-11H-2237-01"