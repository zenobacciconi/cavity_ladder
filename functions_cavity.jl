using ITensors




# light-matter interaction parameter is a global variable
global const var1 = Ref{Float64}(0.0)
function change_g(x::Float64)
    global var1[] = x
    return nothing
end

#some customized operators
function ITensors.op(::OpName"X",
                      ::SiteType"Qudit" ,
                      s::Index)
    d=dim.(s)
    opmat = zeros(d,d)
    for i=1:d-1
        opmat[i,i+1]= sqrt(i)
    end
    opmat +=opmat'
    opmat = opmat/sqrt(2)
  Op=ITensors.itensor(opmat,prime(s),dag(s))
end
function ITensors.op(::OpName"P",
                      ::SiteType"Qudit" ,
                      s::Index)
    d=dim.(s)
    opmat = zeros(ComplexF64,d,d)
    for i=1:d-1
        opmat[i,i+1]= -1im*sqrt(i)
    end
    opmat +=opmat'
    opmat = opmat/sqrt(2)
  Op=ITensors.itensor(opmat,prime(s),dag(s))
end
function ITensors.op(::OpName"X2",
                      ::SiteType"Qudit" ,
                      s::Index)
    d=dim.(s)
    opmat_1 = zeros(ComplexF64,d+1,d+1)
    for i=1:d
        opmat_1[i,i+1]= sqrt(i)
    end
    opmat_1 +=opmat_1'
    opmat_1 =opmat_1*opmat_1
    opmat =zeros(ComplexF64,d,d)
    opmat=opmat_1[1:d,1:d]
    opmat = opmat/2
  Op=ITensors.itensor(opmat,prime(s),dag(s))
end
function ITensors.op(::OpName"P2",
                      ::SiteType"Qudit" ,
                      s::Index)
    d=dim.(s)
    opmat_1 = zeros(ComplexF64,d+1,d+1)
    for i=1:d
        opmat_1[i,i+1]= -1im*sqrt(i)
    end
    opmat_1 +=opmat_1'
    opmat_1 =opmat_1*opmat_1
    opmat =zeros(ComplexF64,d,d)
    opmat=opmat_1[1:d,1:d]
    opmat = opmat/2
  Op=ITensors.itensor(opmat,prime(s),dag(s))
end
function ITensors.op(::OpName"XP",
                      ::SiteType"Qudit" ,
                      s::Index)
    d=dim.(s)
    Xmat = zeros(ComplexF64,d+1,d+1)
    Pmat = zeros(ComplexF64,d+1,d+1)
    for i=1:d
        Xmat[i,i+1]= sqrt(i)
        Pmat[i,i+1]= -1im*sqrt(i)
    end
    Xmat+=Xmat'
    Pmat+=Pmat'
    opmat_1 = Xmat*Pmat
    opmat =zeros(ComplexF64,d,d)
    opmat=opmat_1[1:d,1:d]
    opmat = opmat/2
  Op=ITensors.itensor(opmat,prime(s),dag(s))
end
function ITensors.op(::OpName"PX",
                      ::SiteType"Qudit" ,
                      s::Index)
    d=dim.(s)
    Xmat = zeros(ComplexF64,d+1,d+1)
    Pmat = zeros(ComplexF64,d+1,d+1)
    for i=1:d
        Xmat[i,i+1]= sqrt(i)
        Pmat[i,i+1]= -1im*sqrt(i)
    end
    Xmat+=Xmat'
    Pmat+=Pmat'
    opmat_1 = Pmat*Xmat
    opmat =zeros(ComplexF64,d,d)
    opmat=opmat_1[1:d,1:d]
    opmat = opmat/2
  Op=ITensors.itensor(opmat,prime(s),dag(s))
end

#peierls phase operator and its complex conjugate 
function ITensors.op(::OpName"Expg",
                      ::SiteType"Qudit" ,
                      s::Index)
    d=dim.(s)
    g = var1[] 
    opmat = M_exact(g,d)
  Op=ITensors.itensor(opmat,prime(s),dag(s))
end
function ITensors.op(::OpName"Exp_g",
                      ::SiteType"Qudit" ,
                      s::Index)
    d=dim.(s)
    g = -var1[]
    opmat = M_exact(g,d)
  Op=ITensors.itensor(opmat,prime(s),dag(s))
end

#need an operator with the same name at every site when calling itensors functions
#just construct some that we are not going to use 
function ITensors.op(::OpName"a",
    ::SiteType"Fermion" ,
    s::Index)
opmat = zeros(2,2)
Op=ITensors.itensor(opmat,prime(s),dag(s))
end
function ITensors.op(::OpName"Cdag",
                      ::SiteType"Qudit" ,
                      s::Index)
    d=dim.(s)
    opmat = zeros(d,d)
    for i=1:d-1
        opmat[i,i+1]= sqrt(i)
    end
  Op=ITensors.itensor(opmat,prime(s),dag(s))
end
function ITensors.op(::OpName"C",
                      ::SiteType"Qudit" ,
                      s::Index)
    d=dim.(s)
    opmat = zeros(d,d)
    for i=1:d-1
        opmat[i+1,i]= sqrt(i)
    end
  Op=ITensors.itensor(opmat,prime(s),dag(s))
end


# exact matrix elements of the peierls phase operator exp(i g (a+adag))
function M_exact(g,d)
    L_pol = L_mat(g^2,d)
    M = zeros(d,d).+0im
    i_fact=0.
    j_fact=0.
    sign_g=sign(g)
    for i=1:d
        for j=1:d
            if i>=j
                M[i,j] = (sign_g)^(i-j) *(1im)^(i-j) * exp(-abs(g)^2 /2 + (i-j)*log(abs(g)+10^-10) +j_fact-i_fact) *L_pol[j,i-j+1]
            else
                M[i,j] = (sign_g)^(j-i) *(1im)^(j-i) * exp(-abs(g)^2 /2 + (j-i)*log(abs(g)+10^-10) +i_fact-j_fact) *L_pol[i,j-i+1]
            end
            if abs(M[i,j])<10^-8
                M[i,j]=0.
            end
            j_fact=j_fact+0.5*log(j)
        end
        j_fact=0.
        i_fact=i_fact+0.5*log(i)
    end
    if abs(g)<=10^-10
        M=1. *Matrix(I,d,d)
    end
    return M
end

# function that calculate generalized laguerre polynomials at x=g2 for the two index [alpha,n]
function L_mat(g2,N)
    L_mat=zeros(N,N)
    L_mat[1,:] = ones(N).* 1.
    L_mat[2,:] = (1-g2 .+(0:N-1))
    for j=1:N
        for i=3:N#3:max(N-j+1,3)
            L_mat[i,j]= (2+ ((j-1)-1-g2)/(i-1))* L_mat[i-1,j] - (1+ ((j-1)-1)/(i-1))*L_mat[i-2,j]
        end
    end
    return L_mat
end

# Function that defines the sweeps characteristics, maximum bond dimension and noise parameter.
# Also return a minimum amount of sweeps to reach the required maximum bond dimension
function define_sweeps(chi_fin,N;n_noise=4,n1=4,n2=4,chi0=10)
    #e.g. [20,20,30,30,40,40,60,60,80,80,...]
    #n1=2 , n2=3
    N_req=0
    dchi=10
    trig=true
    chi0_req=chi0
    while chi0_req<chi_fin
	for i=0:n2-1
	    chi0_req+=dchi
	    if chi0_req<chi_fin
		N_req+= n1
	    end 	
	end
	dchi=dchi*2	
    end 
    if N_req>N
	n1=1
	n2=2
	# get [10,20,30,50,70,110,150,230,310,470,630,950,...]
    end
    noise_vec =zeros(N)
    chi_vec = chi_fin*ones(Int,N)
    dchi=10
    n_chi=0
    k=0
    trig=true
    while chi0<chi_fin
        for i=0:n2-1
            if trig
                for j=1:n1
                    chi_vec[k*n1*n2+i*n1+j]=chi0
                end
                if chi0<=chi_fin 
                    trig=true
                end
                chi0+=dchi
                if trig && chi0>chi_fin
                    n_chi = k*n1*n2+i*n1+n1
                    trig=false
                end
            end
        end 
        dchi= dchi*2
        k+=1
    end
    # e.g. [1,1,1e-1,1e-1,1e-2,1e-2,1e-2,0,..]
    # n_chi=7
    # n_noise=3
    # n3=2
    n3=Int(floor(n_chi/n_noise))
    exp=0
    for i=0:n_noise-1
        for j=1:n3
            noise_vec[i*n3+j]= 10. ^-exp
        end
        exp +=1
    end
    k=n3*n_noise
    while k<n_chi
        k+=1
        noise_vec[k]=10. ^(1-exp)
    end
    min_sweep= n_chi
    v1 = vcat(["maxdim"],chi_vec)
    v2 = vcat(["noise"],noise_vec)
    sweep_mat = hcat(v1,v2)
    @show chi_vec,noise_vec
    @show sweep_mat[min_sweep,:]
    return sweep_mat,min_sweep
end


# function that calculate the ground state with dmrg for the ladder+cavity system
function ladder_cavity_gs(g,omega0,t0,t1,t2,sites;flux_ext=0,V=0.,n_sweeps=150,chi=200,use_psi0_prev=false,psi0_prev="initial MPS here",alpha=0.,E_tol=10^-10)
    #g: light-matter coupling
    #t0,t1,t2: hoppings (t2=0 for square ladder)
    #sites: the sites of the MPS must be passed
    #flux_ext: possible external flux
    #V: nn interaction (both vertical and horizontal) 
    #n_sweeps: maximum number of sweeps
    #chi: maximum bond dimension
    #use_psi0_prev,psi0_prev: you can pass an MPS to start with, otherwise the starting one will be a randomly generated
    #apha: symmetry breaking term on the photon of type    alpha*X 
    #E_tol: required accuracy on the energy, after reaching it the dmrg sweeps stops.
    
    

    L2=length(sites)-1 #twice the number of rungs, L2+1 is the actual lenght of the MPS

    #define the sweeps object and the minimum number of sweeps
    #define the initial MPS
    sweeps = Sweeps(n_sweeps)
    if use_psi0_prev  
        setmaxdim!(sweeps,chi)
        psi0=psi0_prev
        min_sweeps=3
    else
        #create an half filled state 
        state = [isodd(n) ? "1" : "0" for n=1:L2+1] 
        #create a random MPS in the half filling sector
        psi0 = randomMPS(sites,state,linkdims=20)
        sweep_mat,min_sweeps=define_sweeps(chi,n_sweeps;n1=2,n2=2,n_noise=n_noise,chi0=20)
        sweeps=Sweeps(n_sweeps,sweep_mat)
    end

    #define the MPO representation of the hamiltonian using the OpSum feature of Itensor
    ampo=OpSum();
    change_g(g/sqrt(L2/2))   #change g as it is a global variable
    phase = exp(1im*flux_ext)    #eventual external flux enters as a phase on the hoppings
    for i=2:L2-1
        if iseven(i)
            #hopping in top leg (cavity dressed)
            ampo+= -t0*phase,"Expg",1,"Cdag",i,"C",i+2 ;
            ampo+= -t0*conj(phase),"Exp_g",1,"Cdag",i+2,"C",i ;
            #hopping top to bottom within one rung
            ampo+= -t1,"Cdag",i,"C",i+1 ;
            ampo+= -t1,"Cdag",i+1,"C",i ;
            #interactions on vertical bond
            ampo+= V,"n",i,"n",i+1 
        else
            #hopping in bottom leg (cavity dressed)
            ampo+= -t0*conj(phase),"Exp_g",1,"Cdag",i,"C",i+2 ;
            ampo+= -t0*phase,"Expg",1,"Cdag",i+2,"C",i ;
            if abs(t2)>10^-8 #only if triangular geometry
                #hopping top to bottom while changing rung
                ampo+= -t2,"Cdag",i,"C",i+1 ;
                ampo+= -t2,"Cdag",i+1,"C",i ;
                #interactions on triangular bond 
                ampo+= V,"n",i,"n",i+1;
            end
        end
        #interaction on horizontal bond
        ampo+= V,"n",i,"n",i+2
    end
    #hopping in the last vertical bond
    ampo+= -t1,"Cdag",L2,"C",L2+1;
    ampo+= -t1,"Cdag",L2+1,"C",L2;
    #interaction in the last vertical bond
    ampo+= V,"n",L2,"n",L2+1;

    #cavity energy
    ampo+= omega0,"n",1;
    if abs(alpha)>10^-10
        #symmetry breaking term
        ampo += alpha,"X",1;
    end

    #define MPO
    H = MPO(ampo,sites);

    #define an observer that stops the dmrg at E_tol and do at least min_sweeps
    obs_dmrg=DMRGObserver(energy_tol=E_tol,minsweeps=min_sweeps)

    #call the 2-site DMRG algorithm of the Itensor library starting from the MPS psi0
    energy,psi = dmrg(H,psi0,sweeps;outputlevel=1,observer=obs_dmrg);

    return energy,psi,obs_dmrg
end



# various functions to compute some observables

function compute_correlation_matrix(psi)
    L2=length(psi)-1
    L=Int(L2/2)

    #built in function of Itensor adjusted for mixed boson-fermion sites
    corr_matsn = my_correlation_matrix(psi,"Cdag","C",sites=2:L2+1)

    corr_mat=zeros(ComplexF64,4,L,L)
    for i=1:L2
        for j=1:L2
            if isodd(i) && isodd(j)
                corr_mat[1,Int((i+1)/2),Int((j+1)/2)]=corr_matsn[i,j]
            elseif isodd(i) && iseven(j)
                corr_mat[3,Int((i+1)/2),Int(j/2)]= corr_matsn[i,j]
            elseif iseven(i) && isodd(j)
                corr_mat[4,Int(i/2),Int((j+1)/2)]= corr_matsn[i,j]
            elseif iseven(i) && iseven(j)
                corr_mat[2,Int(i/2),Int(j/2)]= corr_matsn[i,j]
            end
        end
    end

    return corr_mat
end

#photon entanglement from svd of the MPS
function compute_photon_ent(psi)
    orthogonalize!(psi, 1)
    U,S,V = svd(psi[1], (linkind(psi, 1-1), siteind(psi,1)))
    SvN = 0.0
    for n=1:dim(S, 1)
      p = S[n,n]^2
      SvN -= p * log(p)
    end
    return SvN
end

#photon density matrix encoded as first 20 schimdt states and eigenvalues
function compute_photon_rho(psi;chi0=100)
    orthogonalize!(psi,1)
    U,S,V = svd(psi[1], (linkind(psi, 1-1), siteind(psi,1)))
    chi=dim(S,1)
    N_photon=dim(U,1)
    lambda=zeros(chi0).+10^-16
    vector = zeros(ComplexF64,N_photon,chi0).+10^-16
    for i=1:chi
        lambda[i]=S[i,i]
        for j=1:N_photon
            vector[j,i]=U[j,i]
        end
    end
    return vector[:,1:20],lambda
end


# wigner function given the relevant schmidt states and their eigenvalues
function compute_wigner_function(rho,lam;nx=40,xM=6,pM=6,nP=10,ny=100,yM=20,x0=0)
    x_space=collect(0:nx-1)/nx *xM .-xM/2 .+x0
    y_space=collect(0:ny-1)/ny *yM
    dx=xM/nx
    dy=yM/ny
    nxy= Int(floor(dx/dy+0.1))
    @show nxy
    nY=2*ny+(nx)*nxy
    Y_space=collect(0:nY-1)/nY *(2*ny*dy+(nx)*dx) .-(2*ny*dy+(nx)*dx)/2 .+x0
    p_space=collect(0:nP-1)/nP *pM .-pM/2
    dy=y_space[2]-y_space[1]
    W=zeros(nx,nP)
    @show p_space
    psi_ph=compute_photon_wf(rho,Y_space)
    for ix=1:nx 
        ixy=ny+nxy*ix
        for ip=1:nP
            p=p_space[ip]
            for l=1:10
                W[ix,ip] += sum(real.(exp.(2im*p*y_space) .* conj.(psi_ph[ixy:ixy+ny-1,l]) .*reverse(psi_ph[(ixy-ny+1):(ixy),l])*abs(lam[l])^2))*dy
            end
        end
    end
    return W
end

#real space wavefunction given the expansion in the occupation basis
function compute_photon_wf(photon_rho,x)
    nx=length(x)
    dx=x[2]-x[1]
    N=length(photon_rho[:,1])
    npsi=length(photon_rho[1,:])
    psi=zeros(Complex,nx,npsi)
    exp_gauss=exp.(-0.5*x.^2) .+0*1im
    Herm_0=ones(nx)
    Herm_1=2*copy(x)
    for j=1:npsi
        psi[:,j]+=Herm_0.*exp_gauss/pi^0.25 *photon_rho[1,j]
        psi[:,j]+=sqrt(0.5)*Herm_1.*exp_gauss/pi^0.25 *photon_rho[2,j]
    end


    ifacts=1
    for i=3:N 
        ifacts=ifacts*sqrt(i-1)
        Herm_2=Herm_1.*x *2  -2*(i-2)*Herm_0
        for j=1:npsi 
            psi[:,j]+=Herm_2 .*exp_gauss /2^((i-1)/2)/ifacts/pi^0.25 *photon_rho[i,j]
        end
        Herm_0=copy(Herm_1)
        Herm_1=copy(Herm_2)
    end
    return psi
end

#dressed expectation value < cdag c exp(ig(a+adag)) >
function compute_cdagc_dressed(sites,psi;geometry=1)
    L2=length(psi)-1
    L=Int(L2/2)
    end 
        
    cdagc_dressed=zeros(ComplexF64,2,L)
    vec= [n==1 ? "Expg" : "Id" for n=1:L+1]
    for i=1:L-1
        vec1=copy(vec)
        vec1[2*i]="Cdag"
        vec1[2*i+2]="C"
        J = MPO(sites,vec1)
        cdagc_dressed[1,i]=inner(psi',J,psi)
    end
    vec= [n==1 ? "Exp_g" : "Id" for n=1:L+1]
    for i=1:L-1
        vec1=copy(vec)
        vec1[2*i+1]="Cdag"
        vec1[2*i+3]="C"
        J = MPO(sites,vec1)
        cdagc_dressed[2,i]=inner(psi',J,psi)
    end

    return cdagc_dressed
end


#built in Itensor.correlation_matrix with minor adaptations for the mixed boson-fermion sites
function my_correlation_matrix(psi::MPS, _Op1::AbstractString, _Op2::AbstractString; kwargs...)
    N = length(psi)
    @show 1
    ElT = ComplexF64#promote_itensor_eltype(psi)
    s = siteinds(psi)
  
    Op1 = _Op1 #make copies into which we can insert "F" string operators, and then restore.
    Op2 = _Op2
    onsiteOp = "$Op1 * $Op2"
    fermionic1 = has_fermion_string(Op1, s[2])
    fermionic2 = has_fermion_string(Op2, s[2])
    if fermionic1 != fermionic2
      error(
        "correlation_matrix: Mixed fermionic and bosonic operators are not supported yet."
      )
    end
  
    # Decide if we need to calculate a non-hermitian corr. matrix which is roughly double the work.
    is_cm_hermitian = false #Assume corr-matrix is non-hermitian
    if haskey(kwargs, :ishermitian) #Did the user explicitly request something?
      is_cm_hermitian::Bool = get(kwargs, :ishermitian, false) #Honour users request
    else
      O1 = op(Op1, s, 2)
      O2 = op(Op2, s, 2)
      O1 /= norm(O1)
      O2 /= norm(O2)
      #We need to decide if O1 ∝ O2 or O1 ∝ O2^dagger allowing for some round off errors.
      eps = 1e-10
      is_op_proportional = norm(O1 - O2) < eps
      is_op_hermitian = norm(O1 - dag(swapprime(O2, 0, 1))) < eps
      if is_op_proportional || is_op_hermitian
        is_cm_hermitian = true
      end
      # finally if they are both fermionic and proportional then the corr matrix will
      # be anti symmetric insterad of Hermitian. Handle things like <C_i*C_j>
      # at this point we know fermionic2=fermionic1, but we put them both in the if
      # to clarify the meaning of what we are doing.
      if is_op_proportional && fermionic1 && fermionic2
        is_cm_hermitian = false
      end
    end
  
    if haskey(kwargs, :site_range)
      @warn "The `site_range` keyword arg. to `correlation_matrix` is deprecated: use the keyword `sites` instead"
      sites_ = kwargs[:site_range]
    else
      sites_ = get(kwargs, :sites, 1:N)
    end
    sites = (sites_ isa AbstractRange) ? sites_ : collect(sites_)
  
    start_site = first(sites)
    end_site = last(sites)
  
    psi = copy(psi)
    orthogonalize!(psi, start_site)
    norm2_psi = norm(psi[start_site])^2
  
    # Nb = size of block of correlation matrix
    Nb = length(sites)
  
    C = zeros(ElT, Nb, Nb)
  
    if start_site == 1
      L = ITensor(1.0)
    else
      lind = commonind(psi[start_site], psi[start_site - 1])
      L = delta(dag(lind), lind')
    end
    pL = start_site - 1
  
    for (ni, i) in enumerate(sites[1:(end - 1)])
      while pL < i - 1
        pL += 1
        L = (L * psi[pL]) * dag(prime(psi[pL], "Link"))
      end
  
      Li = L * psi[i]
  
      # Get j == i diagonal correlations
      rind = commonind(psi[i], psi[i + 1])
      C[ni, ni] =
        scalar((Li * op(onsiteOp, s, i)) * prime(dag(psi[i]), not(rind))) / norm2_psi
  
      # Get j > i correlations
      #if !using_auto_fermion() && fermionic2
        Op1 = "$Op1 * F"
      #end
  
      Li12 = (Li * op(Op1, s, i)) * dag(prime(psi[i]))
      pL12 = i
  
      for (n, j) in enumerate(sites[(ni + 1):end])
        nj = ni + n
  
        while pL12 < j - 1
          pL12 += 1
          #if !using_auto_fermion() && fermionic2
         Li12 *= op("F", s[pL12]) * dag(prime(psi[pL12]))
          #else
            #Li12 *= dag(prime(psi[pL12], "Link"))
          #end
          Li12 *= psi[pL12]
        end
  
        lind = commonind(psi[j], Li12)
        Li12 *= psi[j]
  
        val = (Li12 * op(Op2, s, j)) * dag(prime(prime(psi[j], "Site"), lind))
        C[ni, nj] = scalar(val) / norm2_psi
        if is_cm_hermitian
          C[nj, ni] = conj(C[ni, nj])
        end
  
        pL12 += 1
        #if !using_auto_fermion() && fermionic2
         Li12 *= op("F", s[pL12]) * dag(prime(psi[pL12]))
        #else
          #Li12 *= dag(prime(psi[pL12], "Link"))
        #end
        @assert pL12 == j
      end #for j
      Op1 = _Op1 #"Restore Op1 with no Fs"
  
      if !is_cm_hermitian #If isHermitian=false the we must calculate the below diag elements explicitly.
  
        #  Get j < i correlations by swapping the operators
        #if !using_auto_fermion() && fermionic1
          Op2 = "$Op2 * F"
        #end
        Li21 = (Li * op(Op2, s, i)) * dag(prime(psi[i]))
        pL21 = i
        #if !using_auto_fermion() && fermionic1
          Li21 = -Li21 #Required because we swapped fermionic ops, instead of sweeping right to left.
        #end
  
        for (n, j) in enumerate(sites[(ni + 1):end])
          nj = ni + n
  
          while pL21 < j - 1
            pL21 += 1
            #if !using_auto_fermion() && fermionic1
              Li21 *= op("F", s[pL21]) * dag(prime(psi[pL21]))
            #else
             # Li21 *= dag(prime(psi[pL21], "Link"))
            #end
            Li21 *= psi[pL21]
          end
  
          lind = commonind(psi[j], Li21)
          Li21 *= psi[j]
  
          val = (Li21 * op(Op1, s, j)) * dag(prime(prime(psi[j], "Site"), lind))
          C[nj, ni] = scalar(val) / norm2_psi
  
          pL21 += 1
          #if !using_auto_fermion() && fermionic1
           Li21 *= op("F", s[pL21]) * dag(prime(psi[pL21]))
          #else
           # Li21 *= dag(prime(psi[pL21], "Link"))
          #end
          @assert pL21 == j
        end #for j
        Op2 = _Op2 #"Restore Op2 with no Fs"
      end #if is_cm_hermitian
  
      pL += 1
      L = Li * dag(prime(psi[i], "Link"))
    end #for i
  
    # Get last diagonal element of C
    i = end_site
    while pL < i - 1
      pL += 1
      L = (L * psi[pL]) * dag(prime(psi[pL], "Link"))
    end
    lind = commonind(psi[i], psi[i - 1])
    C[Nb, Nb] =
      scalar(L * psi[i] * op(onsiteOp, s, i) * prime(prime(dag(psi[i]), "Site"), lind)) /
      norm2_psi
  
    return C
  end
_op_prod(o1::AbstractString, o2::AbstractString) = "$o1 * $o2"
_op_prod(o1::Matrix{<:Number}, o2::Matrix{<:Number}) = o1 * o
