require("optparse")
options(scipen = 999) # because otherwise it occaisionally prints out scientific notation which screws up stat calculations
  # set up argument parser
  option_list<-list(
    make_option(c("-m", "--map"), type="character", default=NULL,
                help="discoal map file", metavar="<filename>"),
    make_option(c("-p", "--hap"), type="character", default=NULL,
                help="discoal hap file", metavar="<filename>"),
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
    make_option(c("-s", "--start"), type="integer", default=1000000,
			    help="start position of the locus")
  )
  
  opt_parser<-OptionParser(option_list=option_list)
  opt<-parse_args(opt_parser)
  if (is.null(opt$map)){
    print_help(opt_parser)
    stop("You must at least supply the map file (as input file)", call.=FALSE)
    }
  if (is.null(opt$hap)){
    print_help(opt_parser)
    stop("You must at least supply the hap file (as input file)", call.=FALSE)
  }
  if (is.null(opt$outdir)){
    print_help(opt_parser)
    stop("You must also supply the output directory where you want your .hap and .map files", call.=FALSE)
  }


# read map and hap files
print(opt$map)
print(opt$hap)
map<-read.table(opt$map)
hap<-read.table(opt$hap,colClasses="character")
print(hap)
names(map)<-c("chr","ID","genpos","physpos")
map$ID<-rownames(map)
map$physpos<-map$physpos-opt$start

chr<-map[1,1]
    # .map
           # if nsites < 10000, positions won't be an integer

	# restrict to max. window length and center point
	minposition<-opt$center-opt$window/2
	maxposition<-opt$center+opt$window/2
	if(minposition<0 | maxposition>1200000){stop(paste("check center = ",opt$center," and window length = ",opt$window,"\nminposition is ",minposition,", maxposition is ",maxposition))} 

        positions<-map$physpos
	pos_locations<-which(positions>=minposition & positions<=maxposition)
	if(length(pos_locations)<1){stop("no SNPs in this window")}
	phys_positions<-positions[pos_locations]
	gen_positions<-map$genpos[pos_locations]

        locusID<-seq(0,length(phys_positions)-1)
        map_df<-data.frame(chr=rep(chr,length(phys_positions)),locusID,gen_positions,phys_positions)
          # replicate physical position as genetic position
#	print(map_df)

      # write .map file
        basename<-strsplit(opt$map,"[.]")[[1]][1]
	      basename<-strsplit(basename,"/")[[1]][length(strsplit(basename,"/")[[1]])]
        paste(opt$outdir,"/",basename,"_c",as.character(opt$center),".map",sep="")
        write.table(map_df, paste(opt$outdir,"/",basename,"_c",as.character(opt$center),".map",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)

    # .hap
      # make .hap file
#        hap_df<-unite(hap[pos_locations,,drop=FALSE],"snps",everything(),sep="")
        hap_df<-hap[pos_locations,,drop=FALSE]
        print(head(hap))
        print(colnames(hap))
        snps<-apply( hap_df[,c(1:length(colnames(hap)))], 1, paste0, collapse="")
#        hap_df
        hap_isafe<-as.data.frame(hap_df,stringsAsFactors = FALSE)
      #write .hap file
        write.table(snps, paste(opt$outdir,"/",basename,"_c",as.character(opt$center),".hap",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)
      #write .hap file with positions and columns
        write.table(cbind(phys_positions,hap_isafe), paste(opt$outdir,"/",basename,"_c",as.character(opt$center),".hap.iSAFE",sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    
