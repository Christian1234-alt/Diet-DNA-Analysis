/*-------------------------------------------------------------------------*
 *---									---*
 *---		burdickDeconaSummarizer.cpp				---*
 *---									---*
 *---	----	----	----	----	----	----	----	----	---*
 *---									---*
 *---	Version 1a		2026 April 24		Joseph Phillips	---*
 *---									---*
 *-------------------------------------------------------------------------*/

#include	<stdlib.h>
#include	<stdio.h>
#include	<string.h>
#include	<ctype.h>
#include	<unistd.h>
#include	<libgen.h>	// dirname()
#include	<fcntl.h>	// creat()
#include	<errno.h>	// errno
#include	<algorithm>	// std::min
#include	<sys/stat.h>	// stat()
#include	<sys/types.h>	// opendir()
#include	<dirent.h>	// opendir()
#include	<vector>


const char*	FLAG_PREFIXES[]		= { "--",
      					    "-"
					  };

const	size_t	NUM_FLAG_PREFIXES	= sizeof(FLAG_PREFIXES) /
					  sizeof(const char*);

const	char*	HELP_FLAG_STR_ARRAY[]	= { "-?",
					    "-h",
					    "--help"
					  };

const	size_t	NUM_HELP_FLAGS		= sizeof(HELP_FLAG_STR_ARRAY)
					  / sizeof(const char*);

const	int	MAX_LINE		= 4096;

const	char	QUOTE_CHAR		= 34;

const	char	FASTA_ENTRY_BEGIN_CHAR	= '>';

const	char	LOCAL_BLAST_COMMENT_CHAR= '#';

#define		LOCAL_BLAST_FIELD_KEY	" Fields:"

#define		LOCAL_BLAST_HITS_TEMPLATE \
					" %d %s"

#define		CLUSTER_FOLDER_SUFFIX	"_cluster_folder"

#define		CLUSTER_FASTA_FILENAME	"cluster.fasta"

#define		BLAST_OUTPUT_DIR	"blast"

#define		SPECIES_ENTRY_FILEPATH	"blast/remoteSearch.txt"

#define		BLAST_VERSION_PREFIX	"BLAST"

#define		REFERENCE_KEYWORD_PREFIX	\
					"Reference:"

#define		OUTPUT_FILENAME_STR	"results."

#define		LOCAL_PERCENT_ID	"% identity"

#define		LOCAL_SUBJECT_LEN	"subject length"

#define		LOCAL_SUBJECT_TITLE	"subject titles"

#define		LOCAL_GAP_OPENS		"gap opens"

#define		LOCAL_MISMATCHES	"mismatches"

#define		LOCAL_SUBJECT_ID	"subject ids"

#define		QUERY_SUBJECT_LEN	"query length"

const int	NULL_NUM_HITS_VALUE	= -1;

const char*	SPECIES_SPECIMEN_SEPARATOR_ARRAY[]
					= { " voucher ",
					    " isolate ",
					    " strain ",
					    " mitochondrion"
					  };

const	size_t	NUM_SPECIES_SPECIMEN_SEPARATOR
				= sizeof(SPECIES_SPECIMEN_SEPARATOR_ARRAY) /
				  sizeof(const char*);

#define		NIH_LENGTH_PREFIX	"Length="

#define		NIH_IDENTITIES_PREFIX	"Identities"

#define		NIH_GAPS_PREFIX		"Gaps"

const double	MIN_FRACTIONAL_BP_MATCH	= 0.95;

double		minFractionalIdMatch_global
					= MIN_FRACTIONAL_BP_MATCH;

enum	class	BlastOutput :
		int
		{
		  ERROR	= -1,
		  LOCAL,
		  NIH
		};

enum	class	SummarizeOutput :
		int
		{
		  ERROR	= -1,
		  JSON,
		  CSV
		};

const char*	SUMMARIZE_OUTPUT_NAME_ARRAY[]
					= { "json",
					    "csv"
					  };

const	size_t	NUM_SUMMARIZE_OUTPUT_NAME_ARRAY
					= sizeof(SUMMARIZE_OUTPUT_NAME_ARRAY)
					  / sizeof(const char*);

const SummarizeOutput
		DEFAULT_SUMMARIZE_OUTPUT= SummarizeOutput::JSON;

SummarizeOutput	summarizeOutput_global	= DEFAULT_SUMMARIZE_OUTPUT;

char		csvFirstColumnsCArray_global[MAX_LINE];


class	ByteBuffer
{
  //  I.  Member vars:
  //  PURPOSE:  To hold the address of the data to return.
  const char*			dataPtr_;

  //  PURPOSE:  To hold the length of the data to return.
  size_t			numBytes_;

  //  PURPOSE:  To hold the address of the heap-allocated space.
  char*				beginHeapPtr_;

  //  PURPOSE:  To hold the length of the heap-allocated space.
  size_t			numHeapBytes_;

  //  PURPOSE:  To hold the address of where in the heap array to write the
  //		next data.
  char*				heapRunPtr_;

  //  II.  Disallowed auto-generated methods:

protected :
  //  III.  Protected methods:
  //  PURPOSE:  To return the address of just after the heap array.
  const char*	getHeapEndCPtr	()
				const
				{ return(beginHeapPtr_+numHeapBytes_); }

  //  PURPOSE:  To return the remaining bytes left in the heap before it
  //		needs to be resized.
  size_t	getRemainingHeapSize
				()
				const
				{ return(getHeapEndCPtr() - heapRunPtr_); }

public :
  //  IV.  Constructor(s), assignment op(s), factory(s) and destructor:
  //  PURPOSE:  To initialize '*this' to an empty buffer.
  ByteBuffer			() :
				dataPtr_(nullptr),
				numBytes_(0),
				beginHeapPtr_(nullptr),
				numHeapBytes_(0),
				heapRunPtr_(nullptr)
				{ }

  //  PURPOSE:  To release the resources of '*this'.  No parameters.
  //		No return value.
  ~ByteBuffer			()
				{
				  free(beginHeapPtr_);
				}

  //  V.  Accessors:
  //  PURPOSE:  To return a pointer to the beginning of the data.
  //		No parameters.
  const char*	getDataPtr	()
				const
				{
				  return( (dataPtr_==nullptr) ? "" : dataPtr_ );
				}

  //  PURPOSE:  To return the number of bytes of data
  size_t	getNumBytes	()
				const
				{ return(numBytes_); }

  //  PURPOSE:  To return 'true' if '*this' is empty, or 'false' otherwise.
  //		No parameters.
  bool		isEmpty		()
				const
				{ return(getNumBytes() == 0); }

  //  PURPOSE:  To return the string in this stored at offset 'offset'.
  const char*	getCPtr		(size_t		offset
				)
				const
				{
				  return( (offset >= getNumBytes())
				  	  ? ""
					  : (getDataPtr() + offset)
					);
				}

  //  VI.  Mutators:
  //  PURPOSE:  To set '*this' to an empty state.  No parameters.
  //		No return value.
  void		clear		()
				{
				  dataPtr_	= nullptr;
				  numBytes_	= 0;
				  heapRunPtr_	= beginHeapPtr_;
				}

  //  PURPOSE:  To set the byte at index 'offset' to 'b'.  No return value.
  void		setC		(size_t		offset,
		    		 char		b
				)
				{
				  if  ( (offset < getNumBytes())	&&
					(beginHeapPtr_ != nullptr)
				      )
				  {
				    beginHeapPtr_[offset]	= b;
				  }
				}

  //  PURPOSE:  To initialize '*this' to point to the const data pointed to
  //		by 'sourceCPtr' of length 'sourceLen'.  No return value.
  void		initialize	(const char*	sourceCPtr,
				 size_t		sourceLen
				)
				{
				  append(sourceCPtr,sourceLen);
				}

  //  PURPOSE:  To initialize '*this' to point to the const data pointed to
  //		by 'sourceCPtr'.  Length of C-string computed by 'strlen()'.
  //		No return value.
  void		initialize	(const char*	sourceCPtr
				)
				{
				  append(sourceCPtr,strlen(sourceCPtr));
				}

  //  PURPOSE:  To append the C-string pointed to by 'toAppendCPtr' of
  //		length 'length' to the end of '*this'.    Returns offset from
  //		'beginHeapPtr_' where 'toAppendCPtr' will be stored.
  size_t	append		(const char*	toAppendCPtr,
				 size_t		length
				)
  {
    //  I.  Application validity check:

    //  II.  Append 'toAppendCPtr' to the end of '*this':
    size_t	toReturn;

    if  ( (dataPtr_ == nullptr)  ||  (dataPtr_ == beginHeapPtr_) )
    {
      //  II.A.  Handle when there is no const data to copy to heap:
      //  II.A.1.  Allocate space if needed:
      if  ( (numHeapBytes_ == 0)			||
	    (getRemainingHeapSize() < (length + 1))
	  )
      {
	size_t	runOffset	= heapRunPtr_ - beginHeapPtr_;
	size_t	numNeededBytes	= runOffset + length + 1;
	size_t	toAllocate	= (numHeapBytes_ == 0)
	    			  ? MAX_LINE
				  : 2*numHeapBytes_;

	while  (toAllocate < numNeededBytes)
	{
	  toAllocate	*= 2;
	}

	beginHeapPtr_	= (char*)realloc(beginHeapPtr_,toAllocate);

	if  (beginHeapPtr_ == nullptr)
	{
	  fprintf(stderr,"Heap allocation error.\n");
	  exit(EXIT_FAILURE);
	}

	memset(beginHeapPtr_+numHeapBytes_,'\0',toAllocate-numHeapBytes_);
	numHeapBytes_	= toAllocate;
	heapRunPtr_	= beginHeapPtr_ + runOffset;
      }

      //  II.A.2.  Copy 'toAppendCPtr':
      toReturn	 = heapRunPtr_ - beginHeapPtr_;
      memcpy(heapRunPtr_,toAppendCPtr,length);
      heapRunPtr_	+= length;
      *heapRunPtr_	 = '\0';
      numBytes_	 = heapRunPtr_ - beginHeapPtr_;
      dataPtr_	 = beginHeapPtr_;
    }
    else
    {
      //  II.B.  Handle when there *is* const data to copy to heap:
      //  II.B.1.  Remember const data:
      const char*	tempCPtr	= getDataPtr();
      size_t		tempLen		= getNumBytes();

      //  II.B.2.  Reset '*this' so it can hold const or heap data:
      clear();

      //  II.B.3.  Append const data to '*this' in its heap:
      append(tempCPtr,tempLen);

      //  II.B.4.  Concatenate with data originally given to append:
      toReturn	= append(toAppendCPtr,length);
    }

    //  III.  Finished:
    return(toReturn);
  }

  //  PURPOSE:  To append the C-string pointed to by 'toAppendCPtr' to the
  //		end of '*this'.    Returns offset from 'beginHeapPtr_' where
  //		'toAppendCPtr' will be stored.
  size_t	append		(const char*	toAppendCPtr
				)
				{
				  return(append
					    (toAppendCPtr,
					     strlen(toAppendCPtr)
					    )
					);
				}

  //  PURPOSE:  To append the single char 'c' to the end of '*this'.  Returns
  //		offset from 'beginHeapPtr_' where 'toAppendCPtr' will be stored.
  size_t	append		(char	c
				)
				{
				  return(append(&c,1));
				}

  //  PURPOSE:  To append 'rhs' to the end of '*this'.  Returns reference to
  //		'*this'.
  ByteBuffer&	operator+=	(char	rhs
				)
      				{
				  append(rhs);
				  return(*this);
				}

  //  PURPOSE:  To append the C-string pointed to by 'rhsCPtr' to the end of
  //		'*this'.  Returns reference to '*this'.
  ByteBuffer&	operator+=	(const char*	rhsCPtr
				)
      				{
				  append(rhsCPtr);
				  return(*this);
				}

  //  PURPOSE:  To write() the data in '*this' concerning to file descriptor
  //		'fd' with file path 'filePathCPtr'.  Returns num bytes written.
  ssize_t	write		(int		fd,
				 const char*	filePathCPtr
				)
  {
    //  I.  Application validity check:

    //  II.  Attempt to write '*this' to 'fd':
    const size_t  WRITE_BUFFER_LEN	= 0x10000; // 65536
    ssize_t	  ttlNumBytesSent	= 0;	
    ssize_t	  ttlNumBytesToSend	= getNumBytes();

    //  II.A.  Find out about 'fd':

    //  II.B.  Attempt to write '*this' to 'fd':
    while  (ttlNumBytesSent < ttlNumBytesToSend)
    {
      errno			= 0;
      
      size_t	numBytesToSend	= std::min
					(WRITE_BUFFER_LEN,
					 (size_t)
					   (ttlNumBytesToSend-ttlNumBytesSent)
					);
      ssize_t	numBytesJustSent= ::write
					(fd,
					 getDataPtr()+ttlNumBytesSent,
					 numBytesToSend
					);

      if  (numBytesJustSent < 0)
      {
	if  ( (errno == EWOULDBLOCK)  ||  (errno == EAGAIN) )
	{
	  usleep(1000);
	  continue;
	}
	else
	{
	  fprintf(stderr,"write() failed %d (%s)\n",errno,strerror(errno));
	  exit(EXIT_FAILURE);
	}
      }

      ttlNumBytesSent	+= numBytesJustSent;
    }

    //  III.  Finished:
    return(ttlNumBytesSent);
  }

};


class		SpeciesEntry
{
  //  I.  Member vars:
  char*		accessionCPtr_;
  char*		speciesNameCPtr_;
  char*		specimenNameCPtr_;
  char*		geneNameCPtr_;
  int		numMatches_;
  int		numGaps_;
  int		numBasePairs_;
  double	fractionalMatch_;

  //  II.  Disallowed auto-generated methods:

protected :
  //  III.  Protected:

public :
  //  IV.  Constructor(s), assignment op(s), factory(s) and destructor:
  //  PURPOSE:  To initialize '*this' to hold species name 'newSpeciesCPtr',
  //  	gene name 'newGeneCPtr', and factional match 'newFractionalMatch'.
  //	No return value.
  SpeciesEntry			(const char*	newAccessionCPtr,
  				 const char*	newSpeciesCPtr,
  				 const char*	newSpecimenCPtr,
  				 const char*	newGeneCPtr,
				 int		newNumMatches,
				 int		newNumGaps,
				 int		newBasePairs,
				 double		newFractionalMatch
				) :
				accessionCPtr_(strdup(newAccessionCPtr)),
				speciesNameCPtr_(strdup(newSpeciesCPtr)),
				specimenNameCPtr_(strdup(newSpecimenCPtr)),
				geneNameCPtr_(strdup(newGeneCPtr)),
				numMatches_(newNumMatches),
				numGaps_(newNumGaps),
				numBasePairs_(newBasePairs),
				fractionalMatch_(newFractionalMatch)
				{ }

  //  PURPOSE:  To make '*this' a copy of 'sourceRef'.  No return value.
  SpeciesEntry			(const SpeciesEntry&	sourceRef
				) :
				accessionCPtr_
					(strdup(sourceRef.getAccessionCPtr())
					),
				speciesNameCPtr_
					(strdup(sourceRef.getSpeciesNameCPtr())
					),
				specimenNameCPtr_
					(strdup(sourceRef.getSpecimenNameCPtr())
					),
				geneNameCPtr_
					(strdup(sourceRef.getGeneNameCPtr())),
				numMatches_(sourceRef.getNumMatches()),
				numGaps_(sourceRef.getNumGaps()),
				numBasePairs_(sourceRef.getNumBasePairs()),
				fractionalMatch_
					(sourceRef.getFractionalMatch())
				{ }

  //  PURPOSE:  To release the resources of '*this', make '*this' a copy of
  //	'sourceRef', and to return a reference to '*this'.
  SpeciesEntry&	operator=	(const SpeciesEntry&	sourceRef
				)
  {
    if  (this != &sourceRef)
    {
      free(geneNameCPtr_);
      free(specimenNameCPtr_);
      free(speciesNameCPtr_);
      free(accessionCPtr_);

      accessionCPtr_	= strdup(sourceRef.getAccessionCPtr());
      speciesNameCPtr_	= strdup(sourceRef.getSpeciesNameCPtr());
      specimenNameCPtr_	= strdup(sourceRef.getSpecimenNameCPtr());
      geneNameCPtr_	= strdup(sourceRef.getGeneNameCPtr());
      numMatches_	= sourceRef.getNumMatches();
      numGaps_		= sourceRef.getNumGaps();
      numBasePairs_	= sourceRef.getNumBasePairs();
      fractionalMatch_	= sourceRef.getFractionalMatch();
    }

    return(*this);
  }

  //  PURPOSE:  To release the resources of '*this'.
  ~SpeciesEntry			()
				{
				  free(geneNameCPtr_);
				  free(specimenNameCPtr_);
				  free(speciesNameCPtr_);
				  free(accessionCPtr_);
				}

  //  V.  Accessors:
  const char*	getAccessionCPtr()
				const
				{ return(accessionCPtr_); }

  const char*	getSpeciesNameCPtr
				()
				const
				{ return(speciesNameCPtr_); }

  const char*	getSpecimenNameCPtr
				()
				const
				{ return(specimenNameCPtr_); }

  const char*	getGeneNameCPtr	()
				const
				{ return(geneNameCPtr_); }
				
  int		getNumMatches	()
  				const
				{ return(numMatches_); }

  int		getNumGaps	()
				const
				{ return(numGaps_); }

  int		getNumBasePairs	()
				const
				{ return(numBasePairs_); }

  double	getFractionalMatch
				()
				const
				{ return(fractionalMatch_); }

  //  PURPOSE: To output '*this' to 'textCPtr' of length 'textLen' as a JSON
  //	entry.  Returns 'textCPtr'.
  char*		toJson		(char*		textCPtr,
				 size_t		textLen
				)
  {
    snprintf
	(textCPtr,textLen,
	 "\n"
	 "\t{\n"
	 "\t  \"accession\":\"%s\",\n"
	 "\t  \"species\":\"%s\",\n"
	 "\t  \"specimen\":\"%s\",\n"
	 "\t  \"gene\":\"%s\",\n"
	 "\t  \"numMatches\":%d,\n"
	 "\t  \"numGaps\":%d,\n"
	 "\t  \"numBasePairs\":%d,\n"
	 "\t  \"percentId\":%g\n"
	 "\t}",
	 getAccessionCPtr(),
	 getSpeciesNameCPtr(),
	 getSpecimenNameCPtr(),
	 getGeneNameCPtr(),
	 getNumMatches(),
	 getNumGaps(),
	 getNumBasePairs(),
	 getFractionalMatch()
	);
    return(textCPtr);
  }

  //  PURPOSE: To output '*this' to 'textCPtr' of length 'textLen' as a CSV
  //	row.  Returns 'textCPtr'.
  char*		toCsv		(char*		textCPtr,
				 size_t		textLen
				)
  {
    snprintf
	(textCPtr,textLen,
	 "%s\"%s\",\"%s\",\"%s\",\"%s\",%d,%d,%d,%g\n",
	 csvFirstColumnsCArray_global,
	 getAccessionCPtr(),
	 getSpeciesNameCPtr(),
	 getSpecimenNameCPtr(),
	 getGeneNameCPtr(),
	 getNumMatches(),
	 getNumGaps(),
	 getNumBasePairs(),
	 getFractionalMatch()
	);
    return(textCPtr);
  }

};


//  PURPOSE:  To return 'true' if 'pathCPtr' is the path of an existing
//  	ordinary directory or 'false' otherwise.
bool		doesDirectoryExist
				(const char*		pathCPtr
				)
{
  //  I.  Application validity check:

  //  II.  Test for directory existence:
  struct stat	statBuffer;
  int		status		= stat(pathCPtr,&statBuffer);

  //  III.  Return value:
  return( (status==0)  &&  S_ISDIR(statBuffer.st_mode) );
}


FILE*		openDemuxFile	(const char*	demuxFilePathCPtr
				)
{
  FILE*	returnMe	= fopen(demuxFilePathCPtr,"r");

  if  (returnMe == nullptr)
  {
    fprintf
	(stderr,"Cannot open demux file %s, %s.\n",
    	 demuxFilePathCPtr,strerror(errno)
	);
    exit(EXIT_FAILURE);
  }

  return(returnMe);
}


DIR*		openResultsDir	(const char*	resultsDirCPtr
				)
{
  DIR*	returnMe	= opendir(resultsDirCPtr);

  if  (returnMe == nullptr)
  {
    fprintf
	(stderr,"Cannot open results directory %s, %s.\n",
	 resultsDirCPtr,strerror(errno)
	);
    exit(EXIT_FAILURE);
  }

  return(returnMe);
}


int		openOutputFd	(const char*	resultsDirCPtr
				)
{
  char		filePathNameCArray[MAX_LINE];
  size_t	len;
  int		returnMe;
  char*		resultsRealDirCPtr	= realpath(resultsDirCPtr,nullptr);

  strncpy(filePathNameCArray,resultsRealDirCPtr,MAX_LINE);
  filePathNameCArray[MAX_LINE-1]	= '\0';
  dirname(filePathNameCArray);
  free(resultsRealDirCPtr);
  len					= strlen(filePathNameCArray);

  if  ( (len == 0)  ||  (filePathNameCArray[len-1] != '/') )
  {
    filePathNameCArray[len]	= '/';
    len++;
    filePathNameCArray[len]	= '\0';
  }

  strncat(filePathNameCArray+len,OUTPUT_FILENAME_STR,MAX_LINE - len - 1);
  strcat(filePathNameCArray,
	 SUMMARIZE_OUTPUT_NAME_ARRAY[static_cast<int>(summarizeOutput_global)]
	);
  returnMe	= creat(filePathNameCArray,0666);

  if  (returnMe < 0)
  {
    fprintf(stderr,"Cannot open %s, %s.\n",filePathNameCArray,strerror(errno));
    exit(EXIT_FAILURE);
  }

  return(returnMe);
}


bool		doesEndWith	(const char*	textCPtr,
				 const char*	endCPtr
				)
{
  size_t	textLen	= strlen(textCPtr);
  size_t	endLen	= strlen(endCPtr);

  if  (endLen > textLen)
    return(false);

  return(strcmp(textCPtr+(textLen-endLen),endCPtr) == 0);
}


bool		isLineBlank	(const char*	lineCPtr
				)
{
  while  (isspace(*lineCPtr) )
    lineCPtr++;

  return(*lineCPtr == '\0');
}


bool		isAccessionNumber
				(char*		accessionNumber,
				 const char*	sourceCPtr
				)
{
  char*	copyCPtr	= accessionNumber;

  if  ( !isalpha(*sourceCPtr) )
  {
    return(false);
  }

  while  ( isalpha(*sourceCPtr) )
  {
    *copyCPtr++	= *sourceCPtr++;
  }

  if  ( !isdigit(*sourceCPtr) )
  {
    return(false);
  }

  while  ( isdigit(*sourceCPtr) )
  {
    *copyCPtr++	= *sourceCPtr++;
  }

  if  ( isspace(*sourceCPtr)  ||  (*sourceCPtr == '\0') )
  {
    *copyCPtr	= '\0';
    return(true);
  }

  if  (*sourceCPtr != '.')
  {
    return(false);
  }

  *copyCPtr++	= *sourceCPtr++;

  if  ( !isdigit(*sourceCPtr) )
  {
    return(false);
  }

  while  ( isdigit(*sourceCPtr) )
  {
    *copyCPtr++	= *sourceCPtr++;
  }

  if  ( isspace(*sourceCPtr)  ||  (*sourceCPtr == '\0') )
  {
    *copyCPtr	= '\0';
    return(true);
  }

  return(false);
}


bool		didObtainNumSequences
				(int&		numSequences,
				 const char*	dirPathCPtr
				)
{
  char		filePathCArray[MAX_LINE];
  FILE*		filePtr;
  char		line[MAX_LINE];
  bool		isInTheMiddleOfEntry	= false;

  snprintf(filePathCArray,MAX_LINE,"%s/%s",dirPathCPtr,CLUSTER_FASTA_FILENAME);
  numSequences	= 0;
  filePtr	= fopen(filePathCArray,"r");

  if  (filePtr == nullptr)
  {
    return(false);
  }

  while  (fgets(line,MAX_LINE,filePtr) != nullptr)
  {
    if  (isInTheMiddleOfEntry)
    {
      isInTheMiddleOfEntry	= false;
    }
    else
    {
      if  (line[0] == FASTA_ENTRY_BEGIN_CHAR)
      {
	isInTheMiddleOfEntry	= true;
	numSequences++;
      }
    }
  }

  fclose(filePtr);
  return(true);
}


void		removeEndingNewline
				(char*		textCPtr
				)
{
  char*	crCPtr		= strchr(textCPtr,'\r');
  char*	nlCPtr		= strchr(textCPtr,'\n');
  char*	lowerCPtr	= (crCPtr == nullptr)
  			  ? nlCPtr
			  : ( (nlCPtr == nullptr)
			      ? crCPtr
			      : ((crCPtr < nlCPtr) ? crCPtr : nlCPtr)
			    );

  if  (lowerCPtr != nullptr)
  {
    *lowerCPtr	= '\0';
  }
}


char*		getEscapedVersion
				(char*		destCPtr,
				 const char*	sourceCPtr
				)
{
  const char*	LOOK_FOR_CSTR	= "\a\b\f\n\r\t\v\\\"";
  const char*	TRANSLATION_CSTR= "abfnrtv\\\"";
  char*		returnMe	= destCPtr;

  for  ( ;  *sourceCPtr != '\0';  sourceCPtr++)
  {
    const char*	escCPtr	= strchr(LOOK_FOR_CSTR,*sourceCPtr);

    if  (escCPtr == nullptr)
    {
      *destCPtr++	= *sourceCPtr;
    }
    else
    {
      *destCPtr++	= '\\';
      *destCPtr++	= TRANSLATION_CSTR[escCPtr - LOOK_FOR_CSTR];
    }
  }

  *destCPtr	= '\0';
  return(returnMe);
}


BlastOutput	getBlastOutputType
				(FILE*		filePtr
				)
{
  char		line[MAX_LINE];
  BlastOutput	returnMe;

  if  (fgets(line,MAX_LINE,filePtr) == nullptr)
  {
    return(BlastOutput::ERROR);
  }

  returnMe	= (line[0] == LOCAL_BLAST_COMMENT_CHAR)
  		  ? BlastOutput::LOCAL
		  : BlastOutput::NIH;
  rewind(filePtr);
  return(returnMe);
}


void		obtainSpeciesSpecimenAndGene
				(char*			speciesNameCArray,
				 char*			specimenNameCArray,
				 char*			geneNameCArray,
				 const char*		textCPtr
				)
{
  //  I.  Application validity check:

  //  II.  Obtain species, specimen and gene:
  //  II.A.  Initialize names:
  speciesNameCArray[0]	= '\0';
  specimenNameCArray[0]	= '\0';
  geneNameCArray[0]	= '\0';

  //  II.B.  Go past spaces:
  while  ( isspace(*textCPtr) )
  {
    if  (*textCPtr == '\0')
      break;

    textCPtr++;
  }

  if  (*textCPtr == '\0')
  {
    return;
  }

  //  II.C.  Obtain species:
  int		index;
  const char*	foundCPtr	= nullptr;

  for  (index = 0;  index < NUM_SPECIES_SPECIMEN_SEPARATOR;  index++)
  {
    foundCPtr	= strcasestr(textCPtr,SPECIES_SPECIMEN_SEPARATOR_ARRAY[index]);

    if  (foundCPtr != nullptr)
      break;
  }

  if  (foundCPtr != nullptr)
  {
    size_t	len	= foundCPtr - textCPtr;

    strncpy(speciesNameCArray,textCPtr,len);
    speciesNameCArray[len]   = '\0';
    textCPtr		    += len;
    textCPtr		    += strlen(SPECIES_SPECIMEN_SEPARATOR_ARRAY[index]);
  }
  else
  {
    char	tempText[MAX_LINE];
    char	tempText2[MAX_LINE];

    tempText[0]		= '\0';
    tempText2[0]	= '\0';
    sscanf(textCPtr,"%s %s",tempText,tempText2);
    snprintf(speciesNameCArray,MAX_LINE,"%s %s",tempText,tempText2);
    textCPtr		 = strstr(textCPtr,tempText2);
    textCPtr		+= strlen(tempText2);
  }

  //  Get specimen name:
  char*	copyCPtr	= specimenNameCArray;

  if  ( (*textCPtr != ',')  &&  (*textCPtr != '\0') )
  {
    while  ( !isspace(*textCPtr)  &&  (*textCPtr != '\0') )
      *copyCPtr++	= *textCPtr++;

    *copyCPtr	= '\0';

    //  Get gene name:
    copyCPtr	= geneNameCArray;

    while  ( isspace(*textCPtr) )
      textCPtr++;

    while  ((*textCPtr != ',') && (*textCPtr != '\0'))
      *copyCPtr++	= *textCPtr++;

    *copyCPtr	= '\0';
  }

  //  III.  Finished:
}


void		considerAdding	(std::vector<SpeciesEntry>&
							list,
				 const char*		accessionCPtr,
				 const char*		speciesNameCPtr,
				 const char*		specimenNameCPtr,
				 const char*		geneNameCPtr,
				 int   			numMatches,
				 int			numGaps,
				 int			numBasePairs,
				 double			percentIdMatch
				)
{
  size_t	index;
  double	fractionalMatch	= percentIdMatch / 100.0;

  if  (fractionalMatch < minFractionalIdMatch_global)
    return;

  for  (index = 0;  index < list.size();  index++)
  {
    SpeciesEntry&	entry	= list[index];

    if  (strcmp(entry.getSpeciesNameCPtr(),speciesNameCPtr) == 0)
    {
      return;
    }
  }

  if  (index >= list.size())
  {
    list.push_back(SpeciesEntry
			(accessionCPtr,
			 speciesNameCPtr,
			 specimenNameCPtr,
			 geneNameCPtr,
			 numMatches,
			 numGaps,
			 numBasePairs,
			 fractionalMatch
			)
		  );
  }
}


/*---			Local format reading code:			---*/

char*		obtainVersionFromLocal
				(char*		versionCPtr,
				 FILE*		filePtr
				)
{
  char		line[MAX_LINE];
  char*		cPtr;

  *versionCPtr	= '\0';

  do
  {
    if  (fgets(line,MAX_LINE,filePtr) == nullptr)
      break;

    cPtr	= strcasestr(line,BLAST_VERSION_PREFIX);

    if  (cPtr != nullptr)
    {
      strncpy(versionCPtr,cPtr,MAX_LINE);
      removeEndingNewline(versionCPtr);
      break;
    }
  }
  while  (isLineBlank(line));

  return(versionCPtr);
}


bool		didReadLocalFields
				(std::vector<char*>&	fieldsDs,
				 FILE*			filePtr
				)
{
  char	line[MAX_LINE];

  while  (fgets(line,MAX_LINE,filePtr) != nullptr)
  {
    char*	cPtr;

    removeEndingNewline(line);

    if  ( (line[0] == LOCAL_BLAST_COMMENT_CHAR)				&&
    	  ( (cPtr = strcasestr(line+1,LOCAL_BLAST_FIELD_KEY)) != nullptr)
	)
    {
      char*	tokenCPtr;

      for  ( tokenCPtr  = strtok(cPtr + sizeof(LOCAL_BLAST_FIELD_KEY)-1,",");
      	     tokenCPtr != nullptr;
	     tokenCPtr  = strtok(nullptr,",")
	   )
      {
        while  ( isspace(*tokenCPtr) )
	{
	  tokenCPtr++;
	}

	fieldsDs.push_back(strdup(tokenCPtr));
      }

      return(true);
    }
  }

  return(false);
}


bool		didReadLocalNumHits
				(int&		numHits,
				 FILE*		filePtr
				)
{
  char	line[MAX_LINE];
  char	word[MAX_LINE];

  while  (fgets(line,MAX_LINE,filePtr) != nullptr)
  {
    if  ( (line[0] == LOCAL_BLAST_COMMENT_CHAR)				&&
    	  (sscanf(line+1,LOCAL_BLAST_HITS_TEMPLATE,&numHits,word) == 2)	&&
	  ( (strcasecmp(word,"hit")  == 0)  ||
	    (strcasecmp(word,"hits") == 0)
	  )
	)
    {
      return(true);
    }
  }

  return(false);
}


int		getLocalIndex	(std::vector<char*>&	fieldsDs,
				 const char*		fieldNameCPtr
				)
{
  int	index;

  for  (index = 0;  index < (int)fieldsDs.size();  index++)
  {
    if  (strcmp(fieldNameCPtr,fieldsDs[index]) == 0)
      return(index);
  }

  return(-1);
}


void		obtainValues	(std::vector<char*>&		valueDs,
				 char*				lineCPtr
				)
{
  char*	tokenCPtr;

  removeEndingNewline(lineCPtr);

  for  ( tokenCPtr  = strtok(lineCPtr,"\t");
	 tokenCPtr != nullptr;
	 tokenCPtr  = strtok(nullptr,"\t")
       )
  {
    valueDs.push_back(tokenCPtr);
  }

}


std::vector<SpeciesEntry>*
		getSpeciesListFromLocal
				(char*		versionCPtr,
				 char*		referenceCPtr,
				 int&		numHits,
				 FILE*		filePtr
				)
{
  char				line[MAX_LINE];
  std::vector<SpeciesEntry>*	returnMe;
  std::vector<char*>		fieldsDs;

  returnMe			= new std::vector<SpeciesEntry>;

  obtainVersionFromLocal(versionCPtr,filePtr);
  referenceCPtr[0]	= '\0';

  if  ( !didReadLocalNumHits(numHits,filePtr) )
  {
    numHits	= NULL_NUM_HITS_VALUE;
  }

  rewind(filePtr);

  if  ( didReadLocalFields(fieldsDs,filePtr) )
  {
    int		percentIdIndex	= getLocalIndex(fieldsDs,LOCAL_PERCENT_ID);
    int		subjLenIndex	= getLocalIndex(fieldsDs,LOCAL_SUBJECT_LEN);
    int		subjTitleIndex	= getLocalIndex(fieldsDs,LOCAL_SUBJECT_TITLE);
    int		subjIdIndex	= getLocalIndex(fieldsDs,LOCAL_SUBJECT_ID);
    int		gapIndex	= getLocalIndex(fieldsDs,LOCAL_GAP_OPENS);
    int		mismatchesIndex	= getLocalIndex(fieldsDs,LOCAL_MISMATCHES);
    int		queryLenIndex	= getLocalIndex(fieldsDs,QUERY_SUBJECT_LEN);

    if  ( (percentIdIndex >= 0)	&& (subjLenIndex  >= 0) &&
	  (subjTitleIndex >= 0) && (subjIdIndex	  >= 0) &&
	  (gapIndex	  >= 0) && (queryLenIndex >= 0)	&&
	  (mismatchesIndex>= 0)
	)
    {
      while  (fgets(line,MAX_LINE,filePtr) != nullptr)
      {  
	if  (line[0] != LOCAL_BLAST_COMMENT_CHAR)
	{
	  std::vector<char*>	valueDs;
	  char			accessionCArray[MAX_LINE];
	  char			speciesNameCArray[MAX_LINE];
	  char			specimenNameCArray[MAX_LINE];
	  char			geneNameCArray[MAX_LINE];
	  int			numBasePairs	= 0;
	  int			numMatches	= 0;
	  int			numGaps		= 0;
	  double		percentId;
	  char*			c0Ptr;
	  char*			c1Ptr;

	  obtainValues(valueDs,line);

	  if  ( valueDs.size() <=
	      	std::max({percentIdIndex,subjLenIndex,subjTitleIndex,
			  subjIdIndex,gapIndex,mismatchesIndex,queryLenIndex
			 }
			)
	      )
	  {
	    continue;
	  }

	  percentId	= strtod(valueDs[percentIdIndex],nullptr);
	  numGaps	= strtol(valueDs[gapIndex],nullptr,10);
	  numBasePairs	= strtol(valueDs[subjLenIndex],nullptr,10);
	  numMatches	= strtol(valueDs[queryLenIndex],nullptr,10)
	  		  - strtol(valueDs[mismatchesIndex],nullptr,10);

	  if  (isAccessionNumber(accessionCArray,valueDs[subjIdIndex]) )
	  {
	    strncpy(accessionCArray,valueDs[subjIdIndex],MAX_LINE);
	  }
	  else
	  {
	    c0Ptr	= strchr(valueDs[subjIdIndex],'|');

	    if  (c0Ptr != nullptr)
	    {
	      c0Ptr++;

	      c1Ptr	= strchr(c0Ptr,'|');

	      if  (c1Ptr != nullptr)
	      {
		*c1Ptr	= '\0';
	      }

	      if  ( !isAccessionNumber(accessionCArray,c0Ptr) )
	      {
		accessionCArray[0]	= '\0';
	      }
	    }
	    else
	    {
	      accessionCArray[0]	= '\0';
	    }
	  }

	  obtainSpeciesSpecimenAndGene
			(speciesNameCArray,
			 specimenNameCArray,
			 geneNameCArray,
			 valueDs[subjTitleIndex]
			);
	  considerAdding(*returnMe,
			 accessionCArray,
			 speciesNameCArray,
			 specimenNameCArray,
			 geneNameCArray,
			 numMatches,
			 numGaps,
			 numBasePairs,
			 percentId
			);
      	}
      }
    }

    for  (auto& field : fieldsDs)
    {
      free(field);
    }
  }

  return(returnMe);
}


/*---			NIH format reading code:			---*/

char*		obtainVersionFromNih
				(char*		versionCPtr,
				 FILE*		filePtr
				)
{
  char		line[MAX_LINE];  

  *versionCPtr	= '\0';

  do
  {
    if  (fgets(line,MAX_LINE,filePtr) == nullptr)
      break;

    if  (strncasecmp(line,BLAST_VERSION_PREFIX,sizeof(BLAST_VERSION_PREFIX)-1)
	 == 0
	)
    {
      strncpy(versionCPtr,line,MAX_LINE);
      removeEndingNewline(versionCPtr);
      break;
    }
  }
  while  (isLineBlank(line));

  return(versionCPtr);
}


char*		obtainRefFromNih(char*		referenceCPtr,
				 FILE*		filePtr
				)
{
  char		line[MAX_LINE];

  *referenceCPtr	= '\0';

  do
  {
    if  (fgets(line,MAX_LINE,filePtr) == nullptr)
      break;

    if  (strncasecmp
		(line,
		 REFERENCE_KEYWORD_PREFIX,
		 sizeof(REFERENCE_KEYWORD_PREFIX)-1
		)
	 == 0
	)
    {
      const char*	cPtr;

      cPtr  = line + (sizeof(REFERENCE_KEYWORD_PREFIX)-1);

      while  ( isspace(*cPtr) )
        cPtr++;

      strncpy(referenceCPtr,cPtr,MAX_LINE);
      removeEndingNewline(referenceCPtr);

      while  (fgets(line,MAX_LINE,filePtr) != nullptr)
      {
        size_t	len	= strlen(referenceCPtr);

	removeEndingNewline(line);

        if  (isLineBlank(line))
	  break;

	if  (len < (MAX_LINE-3))
	{
	  referenceCPtr[len]	= ' ';
	  len++;
	  referenceCPtr[len]	= '\0';
	}

	strncat(referenceCPtr,line,MAX_LINE-len-1);
	referenceCPtr[MAX_LINE-1]	= '\0';
      }

      break;
    }
  }
  while  (isLineBlank(line));

  return(referenceCPtr);
}


std::vector<SpeciesEntry>*
		getSpeciesListFromNih
				(char*		versionCPtr,
				 char*		referenceCPtr,
				 FILE*		filePtr
				)
{
  char				line[MAX_LINE];
  char				nextLine[MAX_LINE];
  std::vector<SpeciesEntry>*	returnMe;

  returnMe			= new std::vector<SpeciesEntry>;

  obtainVersionFromNih(versionCPtr,filePtr);
  obtainRefFromNih(referenceCPtr,filePtr);

  while  (fgets(line,MAX_LINE,filePtr) != nullptr)
  {
    if  (line[0] == FASTA_ENTRY_BEGIN_CHAR)
    {
      char	accessionCArray[MAX_LINE];
      char	speciesNameCArray[MAX_LINE];
      char	specimenNameCArray[MAX_LINE];
      char	geneNameCArray[MAX_LINE];

      //  Include next line, if is a continuation:
      if  (fgets(nextLine,MAX_LINE,filePtr) == nullptr)
	continue;

      if  (isalpha(nextLine[0]) && (strcasecmp(nextLine,NIH_LENGTH_PREFIX)!=0))
      {
	removeEndingNewline(line);
	strncat(line,nextLine,MAX_LINE-1-strlen(line));
      }

      removeEndingNewline(line);

      //  Get accession number:
      const char*	cPtr = line + 1;

      if  ( !isAccessionNumber(accessionCArray,cPtr) )
      {
	continue;
      }

      obtainSpeciesSpecimenAndGene
			(speciesNameCArray,
			 specimenNameCArray,
			 geneNameCArray,
			 cPtr + strlen(accessionCArray)
			);

      int	numMatches		= 0;
      int	numGaps			= 0;
      int	numBasePairs[2]		= {0,0};
      double	percentIdMatch		= 0.0;

      while  (fgets(line,MAX_LINE,filePtr) != nullptr)
      {
	cPtr	= strcasestr(line,NIH_IDENTITIES_PREFIX);

	if  (cPtr != nullptr)
	{
	  cPtr	+= sizeof(NIH_IDENTITIES_PREFIX)-1;
	  cPtr	 = strchr(cPtr,'=');

	  if  (cPtr != nullptr)
	  {
	    cPtr++;

	    sscanf(cPtr,"%d/%d (%lf%%)",
		   &numMatches,&numBasePairs[0],&percentIdMatch
		  );
	    cPtr	= strcasestr(line,NIH_GAPS_PREFIX);

	    if  (cPtr != nullptr)
	    {
	      cPtr	+= sizeof(NIH_GAPS_PREFIX)-1;
	      cPtr	 = strchr(cPtr,'=');

	      if  (cPtr != nullptr)
	      {
		cPtr++;

		sscanf(cPtr,"%d/%d",&numGaps,&numBasePairs[1]);
	      }
	    }
	  }

	  break;
	}
      }

      considerAdding
		(*returnMe,
		 accessionCArray,
		 speciesNameCArray,
		 specimenNameCArray,
		 geneNameCArray,
		 numMatches,
		 numGaps,
		 numBasePairs[0],
		 percentIdMatch
		);
    }
  }

  return(returnMe);
}


std::vector<SpeciesEntry>*
		getSpeciesList	(char*		versionCPtr,
				 char*		referenceCPtr,
				 int&		numHits,
				 const char*	dirPathCPtr
				)
{
  char				blastDirPathCArray[MAX_LINE];
  DIR*				blastDirPtr;
  std::vector<SpeciesEntry>*	returnMe	= nullptr;

  snprintf(blastDirPathCArray,MAX_LINE,"%s/%s",dirPathCPtr,BLAST_OUTPUT_DIR);
  blastDirPtr	= opendir(blastDirPathCArray);

  if  (blastDirPtr != nullptr)
  {
    struct dirent*	blastEntryPtr;

    for  ( blastEntryPtr  = readdir(blastDirPtr);
    	   blastEntryPtr != nullptr;
	   blastEntryPtr  = readdir(blastDirPtr)
	 )
    {
      if  ( doesEndWith(blastEntryPtr->d_name,".txt") )
      {
        break;
      }
    }

    if  (blastEntryPtr != nullptr)
    {
      char	filePathCArray[MAX_LINE];
      FILE*	filePtr;

      snprintf
	(filePathCArray,MAX_LINE,"%s/%s",
	 blastDirPathCArray,blastEntryPtr->d_name
	);
      filePtr	= fopen(filePathCArray,"r");

      if  (filePtr != nullptr)
      {
	switch  (getBlastOutputType(filePtr) )
	{
	case BlastOutput::LOCAL :
	  returnMe	= getSpeciesListFromLocal
					(versionCPtr,
					 referenceCPtr,
					 numHits,
					 filePtr
					);
	  break;

	case BlastOutput::NIH :
	  returnMe	= getSpeciesListFromNih
					(versionCPtr,
					 referenceCPtr,
					 filePtr
					);
	  numHits	= NULL_NUM_HITS_VALUE;
	  break;

	default :
	  break;
	}

	fclose(filePtr);
      }
    }

    closedir(blastDirPtr);
  }

  return(returnMe);
}


void		handleClusterFolder
				(ByteBuffer&	buffer,
				 const char*	clusterNameCPtr,
				 const char*	sampleNameCPtr,
				 const char*	resultsDirCPtr
				)
{
  char		path[MAX_LINE];
  int		numSequences;
  static
  bool		isFirstCluster	= true;

  snprintf
	(path,MAX_LINE,"%s/%s/%s",
	 resultsDirCPtr,sampleNameCPtr,clusterNameCPtr
	);

  if  ( didObtainNumSequences(numSequences,path) )
  {
    char	text[MAX_LINE];
    char	clusterName[MAX_LINE];
    char	version[MAX_LINE];
    char	reference[MAX_LINE];
    char	referenceEscaped[MAX_LINE];
    int		numHits;
    std::vector<SpeciesEntry>*
		speciesListPtr	= getSpeciesList
					(version,
					 reference,
					 numHits,
					 path
					);
    size_t	nameLen		= strlen(clusterNameCPtr)
				  - (sizeof(CLUSTER_FOLDER_SUFFIX)-1);

    strncpy(clusterName,clusterNameCPtr,nameLen);
    clusterName[nameLen]	= '\0';

    switch  (summarizeOutput_global)
    {
    case SummarizeOutput::JSON :
      snprintf
	(text,MAX_LINE,
	 "%s\n"
	 "  {\n"
	 "    \"sample\":\"%s\",\n"
	 "    \"cluster\":\"%s\",\n"
	 "    \"numSequences\":%d,\n"
	 "    \"blastVersion\":\"%s\",\n"
	 "    \"reference\":\"%s\",\n"
	 "    \"speciesList\":\n"
	 "      [",
	 (isFirstCluster ? "" : ","),sampleNameCPtr,clusterName,numSequences,
	 version,getEscapedVersion(referenceEscaped,reference)
	);
      break;

    case SummarizeOutput::CSV :
      snprintf
	(csvFirstColumnsCArray_global,MAX_LINE,
	 "\"%s\",\"%s\",%d,%d,\"%s\",\"%s\",",
	 sampleNameCPtr,clusterName,numHits,numSequences,
	 version,getEscapedVersion(referenceEscaped,reference)
	);
      text[0]		 = '\0';
      break;
    }
    
    buffer		+= text;

    if  (speciesListPtr != nullptr)
    {
      switch  (summarizeOutput_global)
      {
      case SummarizeOutput::JSON :
	for  (size_t index = 0;  index < speciesListPtr->size();  index++)
	{
	  if  (index > 0)
	  {
	    buffer	+= ',';
	  }

	  buffer	+= (*speciesListPtr)[index].toJson(text,MAX_LINE);
        }
	break;

      case SummarizeOutput::CSV :
	for  (size_t index = 0;  index < speciesListPtr->size();  index++)
	{
	  buffer	+= (*speciesListPtr)[index].toCsv(text,MAX_LINE);
        }
	break;
      }
    }

    switch  (summarizeOutput_global)
    {
    case SummarizeOutput::JSON :
      buffer		+= "\n      ]\n  }";
      break;

    case SummarizeOutput::CSV :
      break;
    }

    isFirstCluster	 = false;
    delete(speciesListPtr);
  }
}


void		handleOneSample	(ByteBuffer&	buffer,
				 const char*	sampleNameCPtr,
				 DIR*  		resultsDirPtr,
				 const char*	resultsDirCPtr
				)
{
  struct dirent*	resultsEntryPtr;
  bool	 		resultsDirNameEndsWithSlash;

  resultsDirNameEndsWithSlash	= (resultsDirCPtr[strlen(resultsDirCPtr)-1]
  				   == '/'
				  );
  rewinddir(resultsDirPtr);

  while  ( (resultsEntryPtr = readdir(resultsDirPtr)) != nullptr )
  {
    if  (strcmp(resultsEntryPtr->d_name,sampleNameCPtr) == 0)
    {
      char		inSampleDirCArray[MAX_LINE];
      DIR*		inSampleDirPtr;
      struct dirent*	inSampleEntryPtr;

      if  (resultsDirNameEndsWithSlash)
	snprintf
	    (inSampleDirCArray,MAX_LINE,"%s%s",
	     resultsDirCPtr,sampleNameCPtr
	    );
      else
	snprintf
	    (inSampleDirCArray,MAX_LINE,"%s/%s",
	     resultsDirCPtr,sampleNameCPtr
	    );

      inSampleDirPtr	= opendir(inSampleDirCArray);

      if  (inSampleDirPtr == nullptr)
      {
	fprintf
	    (stderr,"Cannot open directory %s, %s.\n",
	     inSampleDirCArray,strerror(errno)
	    );
	exit(EXIT_FAILURE);
      }

      while  ( (inSampleEntryPtr = readdir(inSampleDirPtr)) != nullptr )
      {
	if  (doesEndWith(inSampleEntryPtr->d_name,CLUSTER_FOLDER_SUFFIX) )
	{
	  handleClusterFolder
			(buffer,
			 inSampleEntryPtr->d_name,
			 sampleNameCPtr,
			 resultsDirCPtr
			);
	}
      }

      closedir(inSampleDirPtr);
      break;
    }
  }
}


bool		isHelpFlag	(const char*	textCPtr
				)
{
  for  (size_t index = 0;  index < NUM_HELP_FLAGS;  index++)
  {
    if  (strcasecmp(textCPtr,HELP_FLAG_STR_ARRAY[index]) == 0)
      return(true);
  }

  return(false);
}


void		showHelp	(FILE*		outPtr
				)
{
  fprintf
    (outPtr,
     "Usage:\tburdickDeconaSummarizer demuxCsvPath resultsPath [thres] "
     "[%s%s|%s%s]\n"
     "  where:\n"
     "\t'demuxCsvPath' is the path of the CSV file with the sample names,\n"
     "\t    barcodes, and other info on the samples.\n"
     "\t'resultsPath' is the path of the directory with the output dirs on\n"
     "\t    each sample.\n"
     "\t'thres' is the optional threshold of the base-pair match that a gene\n"
     "\t    must have to be reported.  If it is in (1.0 - 100] it is treated\n"
     "\t    as a percent.  If it is in [0.0 - 1.0]"
     " it is treated as a fraction.\n"
     "\t    The default value is %g.\n"
     "\t'%s%s' means output will be in JSON %s\n"
     "\t'%s%s' means output will be in CSV %s\n"
     "\tThis program reads 'demuxCsvPath' to find the names of the samples.\n"
     "\tThen, for each sample in 'resultsPath', it looks for generated\n"
     "\tcluster subdirectories.  Lastly, it writes the file \"%s%s\" or\n"
     "\t\"%s%s\" to the parent directory of 'resultsPath' with its\n"
     "\tfindings.\n",
     FLAG_PREFIXES[0],SUMMARIZE_OUTPUT_NAME_ARRAY[0],
     FLAG_PREFIXES[0],SUMMARIZE_OUTPUT_NAME_ARRAY[1],
     MIN_FRACTIONAL_BP_MATCH,
     FLAG_PREFIXES[0],SUMMARIZE_OUTPUT_NAME_ARRAY[0],
     ( (DEFAULT_SUMMARIZE_OUTPUT == SummarizeOutput::JSON)
       ? "(this is the default)"
       : ""
     ),


  closedir(resultDirPtr);
  fclose(demuxFilePtr);
  return(EXIT_SUCCESS);
}
