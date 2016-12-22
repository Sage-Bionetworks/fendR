
##this is the forest implemetnation of the fendR class


######################################################################
# Create the forestFendR class
#
# This is used to represent a forest-based implementation of the fendR framework

#' An S4 class to represent a Forest-based implementation of the fendR predictive
#' network algorithm
#' @inheritParams fendR
#' @slot forestPath Path to python forest scripts
forestFendR<-setClass("forestFendr",

  #add one more slot for path to python code
  slots=c(
    forestPath="character"
  ),

  prototype=list(forestPath='.'),

  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
  validity=function(object)
  {
    #test to make sure path and correct files exist
  },

  contains="fendR" ##here is where we inherit!

)
