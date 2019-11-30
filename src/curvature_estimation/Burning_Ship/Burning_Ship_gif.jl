using Pkg
pkg"activate ..\\..\\.."

using SymPy
using Plots
using Interpolations
using Images
using Statistics

#=
Choise of function
    BurningShip → f(z) = (|Re(z)| + i|Im(z)|)^2 + c
=#
BurningShip(z, c) = (abs(real(z)) + abs(imag(z))*im)^2 + c

println( BurningShip(Sym("z"), Sym("c")) )
println( BurningShip(Sym("x") + Sym("y")*im, Sym("a") + Sym("b")*im) )

#=
Choice of colour palette
    RGB space
=#
max_iteration = 100
plt_RGBint = plot(legend = false);

#=
# depresión y la ansiedad
pRed    = convert(Vector{UInt8}, [0, 2, 6, 13, 35, 70, 86, 66])
pGreen  = convert(Vector{UInt8}, [17, 27, 45, 56, 85, 104, 79, 63])
pBlue   = convert(Vector{UInt8}, [27, 34, 50, 63, 96, 105, 63, 44])
=#

#=
# frio
pRed    = convert(Vector{UInt8}, [1, 4, 19, 36, 36, 90, 139, 73])
pGreen  = convert(Vector{UInt8}, [21, 27, 41, 59, 79, 132, 98, 56])
pBlue   = convert(Vector{UInt8}, [30, 35, 52, 67, 86, 128, 68, 62])
=#

#=
# inestabilidad emocional y mental
pRed    = convert(Vector{UInt8}, [19, 38, 43, 77, 115, 141, 172, 115])
pGreen  = convert(Vector{UInt8}, [23, 32, 23, 39, 59, 76, 108, 45])
pBlue   = convert(Vector{UInt8}, [24, 32, 16, 26, 34, 44, 60, 35])
=#

#=
# problemas psicológicos
pRed    = convert(Vector{UInt8}, [29, 51, 67, 89, 163, 37, 134, 58])
pGreen  = convert(Vector{UInt8}, [20, 35, 41, 57, 96, 51, 51, 25])
pBlue   = convert(Vector{UInt8}, [15, 20, 16, 18, 7, 26, 21, 16])
=#

pRed    = convert(Vector{UInt8}, [0, 2, 6, 13, 35, 70, 73, 167])
pGreen  = convert(Vector{UInt8}, [17, 27, 45, 56, 85, 104, 92, 185])
pBlue   = convert(Vector{UInt8}, [27, 34, 50, 63, 96, 105, 113, 194])

X   = 0 : (pi/2 + 0.09) / (length(pRed) - 1) : (pi/2 + 0.09)
XN  = 0 : (pi/2 + 0.09) / (max_iteration - 1) : (pi/2 + 0.09)

Red_int = scale(Interpolations.interpolate(pRed, BSpline(Quadratic(Line(OnGrid())))), X)
Green_int   = scale(Interpolations.interpolate(pGreen, BSpline(Quadratic(Line(OnGrid())))), X)
Blue_int    = scale(Interpolations.interpolate(pBlue, BSpline(Quadratic(Line(OnGrid())))), X)

scatter!(plt_RGBint, X, pRed, markercolor = :red);
scatter!(plt_RGBint, X, pGreen, markercolor = :green);
scatter!(plt_RGBint, X, pBlue, markercolor = :blue);

plot!(plt_RGBint, XN, [Red_int(x)   for x = XN], linecolor = :red);
plot!(plt_RGBint, XN, [Green_int(x) for x = XN], linecolor = :green);
plot!(plt_RGBint, XN, [Blue_int(x)  for x = XN], linecolor = :blue);

savefig(plt_RGBint, "results\\interpolation.png")
display(plt_RGBint)

colour_palette = [RGB{Float64}(Red_int(x)/255, Green_int(x)/255, Blue_int(x)/255) for x = XN];
colour_image = Array{RGB{Float64}}(undef, 100, max_iteration);
for i = 1:100
    colour_image[i, :] = colour_palette
end
save(string("results\\colour_palette_", string(max_iteration), "x100.png"), colour_image)
display(plot(colour_image))

#=
Curvature Estimation
=#
function curvature_estimation(points)
    n = length(points)
    if n < 3
        return 0
    end

    i = 0
    j = 1
    k = [abs(atan((points[k] - points[j+=1]) / (points[j] - points[i+=1]))) for k = 3:1:n]

    return mean(k)
end

#=
Generating fractal image
=#
X_intv = [-2.1, 1.4]
Y_intv = [-2.0, 1.0]

#=
X_intv = [-1.75, -0.5]
Y_intv = [-0.5, 0.25]
=#

pixelsX = 240 #240 480 960 1920 3840 7680 15360
pixelsY = 135 #135 270 540 1080 2160 4320 8640

BurningShip_img = zeros(RGB{Float64}, (pixelsY, pixelsX));

Xeval = X_intv[1] : (X_intv[2] - X_intv[1]) / (pixelsX - 1) : X_intv[2]
Yeval = Y_intv[1] : (Y_intv[2] - Y_intv[1]) / (pixelsY - 1) : Y_intv[2]


anim = @animate for max_iteration = 3:33
    println("Iteracion: ", max_iteration)

    i = 1
    for x = Xeval
        j = 1
        for y = Yeval
            z = BigFloat(0)*im
            c = BigFloat(x) + BigFloat(y)*im
            n = 1

            points = zeros(Complex{BigFloat}, max_iteration)
            points[n] = z

            while n < max_iteration
                z = BurningShip(z, c)
                if isnan(z) || isinf(z)
                    break
                end

                n += 1
                points[n] = z
            end

            k = Float64(curvature_estimation(points[1:n]), RoundDown)
            if isnan(k) || isinf(k)
                println("[",i,",",j,"] ->", k)
            elseif k < 0 || k > pi/2 + 0.09
                println("[",i,",",j,"] OR->", k)
            else
                BurningShip_img[j, i] = RGB{Float64}(Red_int(k)/255, Green_int(k)/255, Blue_int(k)/255)
            end

            j += 1
        end
        i += 1
    end

    plot(BurningShip_img)
end

gif(anim, string("results\\Burning_Ship", string(pixelsX), "x", string(pixelsY), "_fps5.gif"), fps = 5)
