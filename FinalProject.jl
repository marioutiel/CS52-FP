using LinearAlgebra
using Plots

function potential_energy(X)
    # Compute the potential energy of the particle configuration
    n = size(X,2)
    energy = 0.0
    for i in 1:n-1
        for j in i+1:n
            r = norm(X[:,i] - X[:,j])
            energy += 1 / r^2
        end
    end
    return energy
end

function gradient(X)
    # Compute the gradient of the potential energy
    k = size(X,1)
    n = size(X,2)
    G = zeros(k, n)
    for i in 1:n-1
        for j in i+1:n
            z = X[:,i] - X[:,j]
            r = norm(z)
            G[:,i] += -2*z / r^4
            G[:,j] += 2*z / r^4
        end
    end
    return G
end

function gradient_check(f, x, g)
    fx = f(x) # computed function 
    gx = g(x) # putative gradient
    
    h = sqrt(eps(fx))
    xi = copy(x)
    gxd = copy(gx) 
    for i=1:length(x)
        xi[i] += h
        gxd[i] = (f(xi) - fx)/h
        xi[i] = x[i] # reset
    end
    absdiff = abs.(gxd .- gx)
    
    return (g=gx, gfd=gxd, maxdiff=maximum(absdiff), normdiff=norm(gxd - gx))
end

function gradient_descent(X::Matrix, alpha::Float64, tol::Float64, max_iter::Int)
    # X = initial point
    # alpha = step size
    # tol = tolerance for convergence
    # max_iter = maximum number of iterations
    iter = 0
    prev_energy = potential_energy(X)
    while iter < max_iter
        G = gradient(X)
        X -= alpha .* G
        X = X ./ sqrt.(sum(X.^2, dims=1)) # re-normalize positions
        
        current_energy = potential_energy(X)
        if abs(current_energy - prev_energy) < tol
            print("\n\nConverged in ", iter, " iterations")
            break
        end
        prev_energy = current_energy
        iter += 1
    end
    
    return X
end

function random_thomson_problem(k::Int, n::Int, alpha::Float64, tol::Float64, max_iter::Int)
    # Initialize particle positions randomly on the unit sphere
    X = randn(k, n)
    X = X ./ sqrt.(sum(X.^2, dims=1))

    # Gradient Descent
    X = gradient_descent(X, alpha, tol, max_iter)
    
    return X
end

function thomson_problem(X::Matrix, alpha::Float64, tol::Float64, max_iter::Int)
    X = X ./ sqrt.(sum(X.^2, dims=1))
    
    X = gradient_descent(X, alpha, tol, max_iter)
    
    return X
end

function plot_sphere(points::Matrix, title::String)
    # Define sphere meshgrid
    b = range(0, 2π, length=50)
    a = range(0, π, length=50)
    x = cos.(b) * sin.(a)'
    y = sin.(b) * sin.(a)'
    z = ones(length(b)) * cos.(a)'

    # Create plot of sphere
    p = plot(
        x, y, z, seriestype=:surface, color=:lightgrey, opacity=0.8,
        xlabel="x", ylabel="y", zlabel="z",
        title=title
    )

    # Add the points
    scatter!(
        points[1, :], points[2, :], points[3, :],
        markersize=5, color=:red, label="Points"
    )

    # Add lines connecting the points
    for i in 1:size(points, 2)
        for j in i+1:size(points, 2)
            plot!(
                [points[1, i], points[1, j]],
                [points[2, i], points[2, j]],
                [points[3, i], points[3, j]],
                color=:black, linewidth=2, linestyle=:solid, label=""
            )
        end
    end

    # Set camera and axis properties for interactivity
    plot!(p, camera=(30, 30), aspect_ratio=:equal)
    plot!(p, xlims=(-1, 1), ylims=(-1, 1), zlims=(-1, 1))
    plot!(p, legend=:bottomleft)

    return p
end

x = [0 1; 1 1; 1 0]
x = x ./ sqrt.(sum(x.^2, dims=1))
check = gradient_check(potential_energy, x, gradient).maxdiff

print("Gradient Check: ", check)

print("\n\n##### RANDOM THOMPSON PROBLEM #####")
X = random_thomson_problem(3,2, 0.01, 1e-10, 10000)
fx = potential_energy(X)
gx = gradient(X)
print("\nEnergy: ", fx, "\nGrad: ", norm(gx), "\nX1: ", X[:,1], "\nX2: ", X[:,2])

print("\n\n##### THOMPSON PROBLEM #####")

X0 = [0 0; 0.5 1; 0.3 -1]
X = thomson_problem(X0, 0.01,1e-10,6000)
fx = potential_energy(X)
gx = gradient(X)
print("\nEnergy: ", fx, "\nGrad: ", norm(gx), "\nX1: ", X[:,1], "\nX2: ", X[:,2])

print("\n\n##### INITIAL POINT ALREADY A SOLUTION#####")

X = [0.0 0.0; 1 -1; 0.0 0.0]
fx = potential_energy(X)
gx = gradient(X)
print("\nEnergy: ", fx, "\nGrad: ", norm(gx), "\nX1: ", X[:,1], "\nX2: ", X[:,2])

print("\n\n##### RANDOM THOMPSON PROBLEM WITH N=3 #####")

X = random_thomson_problem(3,3, 0.01, 1e-10, 10000)
fx = potential_energy(X)
gx = gradient(X)
print("\nEnergy: ", fx, "\nGrad: ", norm(gx))
print("\n\nX1: ",X[:,1], "\nX2: ", X[:,2], "\nX3: ", X[:,3])
p = plot_sphere(X, "3 Particles Distribution on 3D Sphere")
display(p)

print("\n\n##### RANDOM THOMPSON PROBLEM WITH N=4 #####")

X = random_thomson_problem(3,4, 0.01, 1e-10, 10000)
fx = potential_energy(X)
gx = gradient(X)
print("\nEnergy: ", fx, "\nGrad: ", norm(gx))
print("\n\nX1: ",X[:,1], "\nX2: ", X[:,2], "\nX3: ", X[:,3], "\nX4: ", X[:,4])
p = plot_sphere(X, "4 Particles Distribution on 3D Sphere")
display(p)

print("\n\n##### RANDOM THOMPSON PROBLEM WITH N=5 #####")

X = random_thomson_problem(3,5, 0.01, 1e-10, 10000)
fx = potential_energy(X)
gx = gradient(X)
print("\nEnergy: ", fx, "\nGrad: ", norm(gx))
print("\n\nX1: ",X[:,1], "\nX2: ", X[:,2], "\nX3: ", X[:,3] ,"\nX4: ", X[:,4], "\nX5: ", X[:,5])
p = plot_sphere(X, "5 Particles Distribution on 3D Sphere")
display(p)
