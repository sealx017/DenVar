#'@title Compute KDE and JSD
#' 
#' Compute kernel density estimate (KDE) of a marker expression using the function, dens_univ,
#' compute Jensen Shannon Distance (JSD) between two computed KDEs using the function, jensen_shannon_dist,
#' compute the KDE of all the patients of the dataset using the function, Array_KDE, and 
#' compute Jensen Shannon Distance matrix between all the patients using the function, JSD_matrix.
#' 
#' @param x is a list of marker expression values in different cells of a patient
#' @param ngrids is the number of grids used in KDE, default is m = 1024
#' @param px is the KDE of first marker
#' @param py is the KDE of second marker
#' @param Data is the dataset having one column named "SampleID" with the patient IDs and one column with marker expression values
#' @return The function, dens_univ returns KDE and the grid-points as a list, the function, jensen_shannon_dist returns the JSD
#' between two densities the function, Array_KDE returns KDE (and the grid-points) of all the patients in a form of a 3d array,
#' and the function, JSD_matrix returns the JSD distance between all the images in a matrix form.  
#' @export

dens_univ = function(x, ngrids = 1024){
  n = length(x)
  min_coef = 0
  max_coef = 1
  den = matrix(0, nrow = n, ncol = ngrids)
  den_grid = matrix(0, nrow = n, ncol = ngrids)
  for(i in 1:n){
       temp_vec = c(x[[i]][!is.na(x[[i]])])
       if(length(temp_vec)==0){
         den[i,] = rep(0,ngrids)
         }
      else{
      if(range(temp_vec)[1]==range(temp_vec)[2]){
        den[i,1] = 1
      }
      else if(length(temp_vec)<5){
        s = density(temp_vec, from = min_coef, to = max_coef,
        n = ngrids, bw = mean(range(temp_vec)))
        den[i,] = s$y
        den_grid[i,] = s$x
      }
      else{
        s = density(temp_vec, from = min_coef, to = max_coef,
         n = ngrids, bw='nrd0')
        den[i,] = s$y
        den_grid[i,] = s$x
      }}}
  return(list(den, den_grid))
}

jensen_shannon_dist = function(px,py){
  px = px/sum(px); py = py/sum(py); 
  px[which(px < .Machine$double.xmin)] <- .Machine$double.xmin
  py[which(py < .Machine$double.xmin)] <- .Machine$double.xmin
  pmean = 1/2*(px + py)
  JSD = sqrt(KLD(px, pmean)$sum.KLD.px.py + KLD(py, pmean)$sum.KLD.px.py)
  return(JSD)
} 


Array_KDE = function(Data, ngrids = 1024){
  sel_images = unique(Data$SampleID)
  Array_dens = array(NA, dim = c(length(sel_images), 2, ngrids))
  x = NULL
  which_images = 1
  len_cov = 1
  for(images in sel_images){
   data_image = Data[Data$SampleID==images,]
  for(i in 1:len_cov)
  {
    x[[i]] = na.omit(as.matrix(data_image[,2]))
  }
    den = dens_univ(x, ngrids = ngrids)
    Array_dens[which_images,1,] = den[[1]]
    Array_dens[which_images,2,] = den[[2]]
    which_images = which_images + 1
  }
  return(Array_dens)
  }
  

JSD_matrix = function(Array_dens){
  num_images = dim(Array_dens[,1,])[1]
  rpd_mat = matrix(0, num_images, num_images)
  for( i in 1:num_images){
   for(j in 1:num_images){
    if(j > i){
      rpd_mat[i,j] <- jensen_shannon_dist(Array_dens[i,1,], Array_dens[j,1,])
    }
  }
  }
  rpd_mat <- rpd_mat+t(rpd_mat)
  colnames(rpd_mat) <- rownames(rpd_mat) <- sel_images
  check_na_row = which(is.na(rpd_mat[,1])==T)
  if(length(check_na_row)>0){
   rpd_mat <- rpd_mat[-check_na_row,-check_na_row]}
  rpd_mat2 <- rpd_mat
  rpd_mat2[rpd_mat2>quantile(rpd_mat,0.95)] <- quantile(rpd_mat,0.95)
  JSD_matrix <- rpd_mat2/max(rpd_mat2)
  return(JSD_matrix)
 }