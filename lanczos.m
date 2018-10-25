# Working implementation of Lanczos algorithm.
# Seems to be very good at predicting highest eigenvalue. Moderately good at predicting the lowest eigenvalue. Everything in between is a crapshoot.
# But this is what Lanczos is supposed to do anyway.
# Written to make sure I understood what was going on with the Lanczos algorithm, so that I could do the implementation in C.

L = 4
n = 2^L
m = 6 # Dimension of Krylow space.
A = rand(n,n) # Random matrix.
A = (A+A')/2 # Make matrix symmetric.
T = zeros(m,m) # Krylow matrix.
v0 = rand(n,1) # Random vector.
v0 = v0/norm(v0) # Normalize.
w = A*v0
T(1,1) = v0'*w
f = w - T(1,1)*v0

# Actual algorithm.

for j = 0:m-2
    b = norm(f) # Get b_j.
    T(j+1,j+2) = b # Put b_j in matrix.
    T(j+2,j+1) = T(j+1,j+2) # Symmetrize.
    v = f/b # Get v_{j+1}.
    w = A*v-b*v0 
    v0 = v # Update v0.
    a = v'*w # Get alpha_{j+1}.
    T(j+2,j+2) = a # Put in matrix.
    f = w - a*v 
    
end

eig(A)
eig(T)
