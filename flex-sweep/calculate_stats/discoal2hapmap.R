## Script to convert discoal (ms) ouput to .hap + .map output

###### consider adapting this to accomodate wildcard expansions or lists of files

#require("gap")

source("calculate_stats/read.ms.output.R")
require("optparse")
options(scipen = 999) # because otherwise it occaisionally prints out scientific notation which screws up stat calculations

  # set up argument parser
  option_list<-list(
    make_option(c("-f", "--file"), type="character", default=NULL,
                help="discoal output file", metavar="<filename>"),
    make_option(c("-o", "--outdir"), type="character", default=NULL,
		help="full path to output directory where you want the .hap and .map files",
		metavar="/path/to/outdir"),
    make_option(c("-l", "--locuslength"), type="integer",
		            help="how many bp in the locus", metavar="# bp", 
			    default=1200000),
    make_option(c("-c", "--center"), type="integer", default=600000,
			    help="where is the center point to use?"),
    make_option(c("-w", "--window"), type="integer", default=1000000,
			    help="how long is the longest window?"),
    make_option(c("-r", "--recombination_map"), type="character", default=NULL,
                help="full path to the recombination map")
  )
  
  opt_parser<-OptionParser(option_list=option_list)
  opt<-parse_args(opt_parser)
  
  if (is.null(opt$file)){
    print_help(opt_parser)
    stop("You must at least supply the discoal output file (as input file)", call.=FALSE)
  }
  if (is.null(opt$outdir)){
    print_help(opt_parser)
    stop("You must also supply the output directory where you want your .hap and .map files", call.=FALSE)
  }

  # read discoal output
    discoal<-read.ms.output(opt$file, is.file=TRUE)
#    discoal
    # assume there is just one replicate run per file
      chr<-1
    # get number of sites (since discoal reports position on a [0,1] scale)
      nsites<-opt$locuslength

    # .map
      # make .map file
         positions<-round(discoal$positions[[1]]*nsites) # this rounds the positions to nearest integer
           # since positions from discoal are in interval [0,1] with five(!) decimal places
           # if nsites < 10000, positions won't be an integer
         duplicates<-duplicated(positions)
         sum(duplicates)
          # duplicates, since rounding above as well as rounding by discoal (when window is large) could make two or more successive SNPs have same position
         while(sum(duplicates)>=1){ # while there are any duplicates left, increment all duplicated positions by 1
             positions[duplicates]<-positions[duplicates]+1
             #positions
             duplicates<-duplicated(positions)
             sum(duplicates)
         }
        positions<-positions+1 # want 1 indexed not 0

            # restrict to max. window length and center point
        minposition<-opt$center-opt$window/2
        maxposition<-opt$center+opt$window/2
	if(minposition<0 | maxposition>1200000){stop(paste("check center = ",opt$center," and window length = ",opt$window,"\nminposition is ",minposition,", maxposition is ",maxposition))} 

        pos_locations<-which(positions>=minposition & positions<=maxposition)
	if(length(pos_locations)<1){
            print("no SNPs in this window")
            map_df<-data.frame(rep(0,4))               
        } else {
            positions<-positions[pos_locations]
	    if (is.null(opt$recombination_map)){
                genpositions<-positions
            } else {
                recomb_map<-read.delim(opt$recombination_map)
                recomb<-data.frame(pos=c(1,recomb_map[,1]),rate=c(recomb_map[1,2],recomb_map[,2]))
                pos<-c()
                rate<-c()

                for (i in 1:(length(recomb[,1])-1)){
#                    if (recomb[i,2]==recomb[i+1,2]){
                        N<-recomb[i+1,1]-recomb[i,1]
                        pos<-c(pos,seq(recomb[i,1],recomb[i+1,1]-1))
                        rate<-c(rate,rep(recomb[i,2],N))
 #                   } else {
  #                      N<-recomb[i+1,1]-recomb[i,1]
   #                     pos<-c(pos,seq(recomb[i,1],recomb[i+1,1]-1))
    #                    rate<-c(rate,rep(recomb[i,2],N))
     #               }
                }        
                i<-i+1
                pos<-c(pos,recomb[i,1])
                rate<-c(rate,recomb[i,2])
                N<-length(pos)
                diff<-c(0,1e-6*rate[-N]*(pos[-1]-pos[-N]))
                all_genpos<-cumsum(diff) + 1 # want 1 indexed
                genpositions<-all_genpos[positions]
            }

        locusID<-seq(0,length(positions)-1)
        map_df<-data.frame(rep(chr,length(positions)),locusID,genpositions,positions)
          # replicate physical position as genetic position
#	print(map_df)
        }
      # write .map file
        basename<-strsplit(opt$file,"[.]")[[1]][1]
	      basename<-strsplit(basename,"/")[[1]][length(strsplit(basename,"/")[[1]])]
        write.table(map_df, paste(opt$outdir,"/",basename,"_c",as.character(opt$center),".map",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)

    # .hap
      # make .hap file
        haplotypes<-apply(discoal$gametes[[1]], 1, paste, collapse="")
        hap_df<-data.frame(haplotypes)
		hap_df<-hap_df[pos_locations,,drop=FALSE]
        hap_isafe<-sapply(haplotypes, strsplit, split="")
		hap_isafe<-t(as.data.frame(hap_isafe,stringsAsFactors = FALSE))[pos_locations,,drop=FALSE]
      #write .hap file
        write.table(hap_df, paste(opt$outdir,"/",basename,"_c",as.character(opt$center),".hap",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)
      #write .hap file with positions and columns
        write.table(cbind(positions,hap_isafe), paste(opt$outdir,"/",basename,"_c",as.character(opt$center),".hap.iSAFE",sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
      #write .hap file with transposed columns
#        write.table(t(hap_isafe), paste(opt$outdir,"/",basename,"_c",as.character(opt$center),".hap.selscan",sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    
