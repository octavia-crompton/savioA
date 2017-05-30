import numpy as np

def make_phi(folder = 'test', fname = 'vanG', dz = 0.2, zmin = 0, zmax = 20, z = 0,
              zs = 4, infL = 2, ksatB = 0.19/3600. , ksatV = 3.38/3600. ):
        # Define van Genuchten parameters
        # veg areas / below the seal layer
        vegAlpha   = 0.0096  
        vegTheta_S = 0.472  
        vegTheta_R = 0.0378  
        vegLambdA  = 0.318 
        vegn       = vegLambdA + 1
        vegm       = vegLambdA/vegn
        ksatV    = ksatV 

        # Seal layer
        sealAlpha   = 0.0078  
        sealTheta_S = 0.450  
        sealTheta_R = 0.0394  
        sealLambdA  = 0.263  
        sealn       = sealLambdA + 1
        sealm       = sealLambdA/sealn
        ksatB    = ksatB
    
        # veg
        nz = z.shape[0]
        alpha = np.ones(nz)*vegAlpha
        theta_S = np.ones(nz)*vegTheta_S
        theta_R = np.ones(nz)*vegTheta_R
        lambdA = np.ones(nz)*vegLambdA
        ksat = np.ones(nz)*ksatV
        
        n = lambdA + 1
        m = lambdA/n
        
        phi_veg = {'alpha': alpha.copy(),
               'theta_R': theta_R.copy(),
               'theta_S': theta_S.copy(),
               'theta_S': theta_S.copy(),
               'lambdA': lambdA.copy(),
               'n': n.copy(),
               'm': m.copy(),
               'ksat': ksat.copy(),
              }
     
        # seal 
        si = np.where(z == z[-1] - zs)[0][0]  
        nz = z.shape[0]

        alpha[si:] = sealAlpha
        theta_S[si:] = sealTheta_S
        theta_R[si:] = sealTheta_R
        lambdA[si:] = sealLambdA
        ksat[si:] = ksatB

        n = lambdA + 1
        m = lambdA/n
        
        phi_bare = {'alpha': alpha,
           'theta_R': theta_R,
           'theta_S': theta_S,
           'theta_S': theta_S,
           'lambdA': lambdA,
           'n': n, 
           'm': m, 
           'ksat': ksat,
          }
        
        return phi_veg, phi_bare
        
        
def write_phi( phi_veg ,phi_bare , hinit,folder = 'test', fname = 'vanG' , 
              infL = 2, dz = 0.2, nz = 101, zmax = 20):
        f = open('{0}/input/{1}.dat'.format(folder, fname), 'w')
        f.write('infL \n')  
        f.write('{0:<13}  \n'.format(infL))
        f.write('dz zmax \n')  
        f.write('{0:<13} {1:<13}  \n'.format(dz, zmax))
        f.write('nz \n')  
        f.write('{0:<13}  \n'.format(nz))        
        
        for ndum in range(nz):
            f.write('{0} {1} {2} {3} {4} {5} \n'.format(phi_veg['alpha'][ndum],
                            phi_veg['theta_S'][ndum], phi_veg['theta_R'][ndum], 
                            phi_veg['lambdA'][ndum], phi_veg['ksat'][ndum], hinit[ndum]))
        
       
        f.write(' bare soil params \n'.format(nz))
        for ndum in range(nz):
            f.write('{0} {1} {2} {3} {4} {5} \n'.format(phi_bare['alpha'][ndum],
                            phi_bare['theta_S'][ndum], phi_bare['theta_R'][ndum], 
                            phi_bare['lambdA'][ndum], phi_bare['ksat'][ndum], hinit[ndum]))
        
        f.close()