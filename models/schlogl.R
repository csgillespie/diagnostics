library(issb)
h = function(x, pars) {
  hazs = numeric(length(pars))
  hazs[1] = pars[1]*x[1]*(x[1]-1)*x[2]/2
  hazs[2] = pars[2]*x[1]*(x[1]-1)*(x[1]-2)/6
  hazs[3] = pars[3]*x[3]
  hazs[4] = pars[4]*x[1]
  return(hazs)
}

#The stoichiometric matrix
## 3 species four reactions
## species 2 & 3 constant ------
smat = matrix(0,nrow=3,ncol=4)
smat[1,1] = 1;
smat[2,1] = -1;
smat[1,2] = -1
smat[2,2] = 1

smat[1,3] = 1
smat[3,3] = -1
smat[1,4] = -1
smat[3,4] = 1
rownames(smat) = c("X", "A", "B")
smat

## The Jacobian
## 4 rows/reactions, 3 columns/species---------
f = get_f = function(x, pars)
{
  fmat = matrix(0, nrow=length(pars), ncol=length(x))
  fmat[1,1] = pars[1]*x[2]/2*(2*x[1]-1)
  fmat[1,2] = pars[1]*x[1]*(x[1]-1)/2
  
  fmat[2,1] = pars[2]*(3*x[1]^2 -6*x[1] + 2)/6
  fmat[3,3] = pars[3]
  fmat[4,1] = pars[4]
  fmat
}

#Build the model
initial = c(250, 1e5, 2e5)
pars = c(3e-7, 1e-4, 1e-3, 3.5)
model = create_model(smat, h, initial, pars, f)
model$get_initial()
model$get_stoic()
