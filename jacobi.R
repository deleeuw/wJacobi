
# compile and link jacobi.c
# load jacobi.dll
# source jacobi function below

jacobi <- function( a, numcomp = ncol( a ), symmetric = FALSE, only.values = FALSE )
# truncated eigenvalue decomposition using jacobi, mimicing eigen()
# function is faster with symmetric matrix, symmetric = TRUE, and numcomp < ncol( a )
{
  a <- unname( as.matrix( a ) )
  n <- nrow( a )
  if ( !n ) stop( "jacobi: 0 x 0 matrix" )
  if ( n != ncol( a ) ) stop( "jacobi: non-square matrix in 'jacobi'" )
  n <- as.integer( n )
  if ( is.na( n ) ) stop( "jacobi: invalid nrow( a )" )
  if ( !all( is.finite( a ) ) ) stop( "jacobi: infinite or missing values in 'a'" )
  if ( missing( symmetric ) ) symmetric <- isSymmetric.matrix( a )
  if ( !symmetric ) stop( "jacobi: non-symmetric matrix in 'jacobi'" ) 
  if ( numcomp < 1 || numcomp > n ) stop( "jacobi: invalid number of components requested" )
  vecs <- rep( 0, n * n )
  vals <- rep( 0, n )
  r <- ( .C( "Cjacobi", n=as.integer(n), a=as.double(a), vecs=as.double(vecs), vals=as.double(vals), numcomp=as.integer(numcomp) ) )
  if ( r$numcomp == 0 ) stop( "jacobi: invalide eigenvalue decomposition" )
  if ( only.values ) return( list( values = r$vals[1:numcomp], vectors = NULL ) )
  else structure( class = "eigen", list( values = r$vals[1:numcomp], vectors = matrix( r$vecs, n, n )[,1:numcomp] ) )
} # jacobi
