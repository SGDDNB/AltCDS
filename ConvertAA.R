#load the libraries required
suppressMessages(library(readr))
suppressMessages(library(Biostrings))
suppressMessages(library(seqinr))

#this will swap and Ts for a U since we are assuming to be working in RNA space
swapTtoU <- function(sequence){
  sequence <- gsub("T", "U", sequence) 
  return(sequence)
}

#this will take a codon and get the corresponding AA from the supplied codon table
convertToAA <- function(codon,table){
  lookup_AA <- table[which(table$codon == codon),]$AA
  return(lookup_AA)
}

#this will take a codon and convert it the most common alternative codon with probability = prob_of_change
convertToAlternativeCodons <- function(codon,table,prob_of_change){
  lookup_AA <- table[which(table$codon == codon),]$AA
  if(runif(1) > (1-prob_of_change)){
    alternative_codon <- getMostFrequentAlternative(codon,lookup_AA,table)
  }else{
    alternative_codon <- codon
  }
  
  return(alternative_codon)
}

#this will take the codon and the AA to the most common alternative
getMostFrequentAlternative <- function(codon,AA,table){
  possibleCodons <- table[which(table$AA == AA),]
  possibleCodons <- possibleCodons[order(-possibleCodons$fraction),]
  number_of_options = dim(possibleCodons)[1]
  if(number_of_options == 1){
    return(codon)
  }
  
  firstAlternative = possibleCodons[1,1]
  secondAlternative = possibleCodons[2,1]
  if(firstAlternative == codon){
    return(secondAlternative)
  }else{
    return(firstAlternative)
  }
}

#this function takes the original sequence and converts it to codons and then gets the the alternative sequence
getAA <- function(sequence, table, prob_of_change){
  AAseq <- list();
  OriginalCodons <- list();
  AlternativeCodons <- list();
  converted_sequence <- swapTtoU(sequence)
  codons = substring(converted_sequence, seq(1, nchar(converted_sequence)-1, 3), seq(3, nchar(converted_sequence), 3))
  number_of_codons = length(codons)
  for (i in 1:number_of_codons){
    AA <- convertToAA(codons[i],table)
    OriginalCodons[i] <- codons[i]
    AlternativeCodons[i] <- convertToAlternativeCodons(codons[i],table,prob_of_change)
    AAseq[i] <- AA
  }
  return(list(AAseq,OriginalCodons,AlternativeCodons))
}


#this function takes a sequence and processes it to get the alternative transcript, it keeps repeating this up to 1000 times if the resulting string contains pallendromes
process_sequence <- function(sequence,CodonUsage,prob){
loops <- 1;
cont = 1;
number_of_palindromes <- 999
min_number_of_palindromes <- 999;
FinalConverted = '';
while(cont == 1){
  NewSequenceInfo = getAA(sequence,CodonUsage,prob)
  AAsequence = paste(NewSequenceInfo[[1]],collapse='')
  Original = paste(NewSequenceInfo[[2]],collapse='')
  Converted = paste(NewSequenceInfo[[3]],collapse='')
  convertedBString <- BString(x=Converted, start=1, nchar=NA)
  palInConverted <- findPalindromes(convertedBString, min.armlength=10,max.looplength=1, min.looplength=0, max.mismatch=0)
  number_of_palindromes <- length(palInConverted)
  if(number_of_palindromes < min_number_of_palindromes){
    min_number_of_palindromes = number_of_palindromes
    FinalConverted = Converted
    if(min_number_of_palindromes == 0){
      cont = 0;
    }
  }
  loops <- loops + 1
  if(loops > 1000){
    cont = 0;
  }
}
return(list(Original,FinalConverted,AAsequence,number_of_palindromes))
}

#this takes a fasta file location and executes the script on each entry
processFile = function(filepath) {
  fasta_data <- read.fasta(file = filepath, as.string = TRUE)
  number_of_seqs = length(fasta_data)
  i = 1
  results = list()
  while ( i <= number_of_seqs ) {
    results[[i]] <- append(process_sequence(toupper(fasta_data[[i]][1]),CodonUsage,0),list(fasta_data[[i]][1],attr(fasta_data[[i]],'name')))
    i = i + 1
  }
  return(results)
}

#this takes the results and prints them
printResults = function(results,filepath) {
  length_of_results = length(results)
  output_file = paste0(filepath,".mimic.fasta")
  if (file.exists(output_file)){
    file.remove(output_file)
  }
  for (seq in 1:length_of_results){
    results_strings = results[[seq]]
    write(paste0(">Alternative (",results_strings[[4]]," pallindromes) for ",results_strings[[6]]),file=output_file,append=TRUE)
    write(results_strings[[2]],file=output_file,append=TRUE)
  }
  return(output_file)
}

#this is the main code of the script - it takes the arguments and executes
#################################CONVERTAA.R##############################
args = commandArgs(trailingOnly=TRUE)
filepath = 'exampleSeqs.txt'
codon_file = "CodonUsage.csv"
if (length(args)==0) {
  stop("Usage: Rscript ConvertAA.R [input.fasta] [OPTIONAL codon.csv, format UUU,F,0.46,17.6,-714298]", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  filepath = args[1]
} else if (length(args)==2){
  filepath = args[1]
  codon_file = args[2]
}
CodonUsage <- read_csv(codon_file, col_names = FALSE,col_types = cols())
colnames(CodonUsage) <- c('codon','AA','fraction','frequency','number')
results <- processFile(filepath)
outfile_path <- printResults(results,filepath)
write(paste("Output be found in: ",outfile_path), stderr())
