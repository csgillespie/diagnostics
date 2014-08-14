library(issb)
h = function(x, pars) {
  hazs = numeric(length(pars))
  hazs[1] = pars[1]
  hazs[2] = pars[2]
  hazs[3] = pars[3]*x[1]
  hazs[4] = pars[4]*x[2]
  hazs[5] = pars[5]*x[1]*x[2]
  return(hazs)
}


#The stoichiometric matrix
smat = matrix(0,nrow=2,ncol=5)
smat[1,1] = 1
smat[2,2] = 1
smat[1,3] = -1
smat[2,4] = -1
smat[1,5] = -1
smat[2,5] = 1
rownames(smat) = c("Y1", "Y2")
##The Jacobian
f = get_f = function(x, pars)
{
  fmat = matrix(0, nrow=5, ncol=2)
  fmat[1,1] = 0
  fmat[2,2] = 0
  fmat[3,1] = pars[3]
  fmat[4,2] = pars[4]
  fmat[5,1] = pars[5]*x[2]
  fmat[5,2] = pars[5]*x[1]
  fmat
}

#Build the model
initial = c(0, 0)
theta = c(10,1,0,0.7,0.008)
model = create_model(smat, h, initial, theta, f)
model$get_initial()
model$get_stoic()
model$get_haz(c(100, 100))
