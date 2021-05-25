"""
Defines the power model of the Huwawei P30 Pro display. The R, G and B arguments
are the raw red, green and blue intensities (between 0 and 255). For the raw
fit values please refer to the power Jupyter Notebook.
"""
# TODO double check the values and finish documentation
function power(R, G, B)
    p_r = 9.33773772e-6.*R.^3 .- 1.77241480e-3.*R.^2 .+ 3.63591812e-1.*R .+ 0.40109990900690207
    p_g = 3.31277854e-6.*G.^3 .+ 9.59555435e-4.*G.^2 .+ 2.57421833e-2.*G .+ 2.9828252112477784
    p_b = 4.19953259e-6.*B.^3 .+ 9.91716987e-4.*B.^2 .+ 1.99992989e-1.*B .+ 7.100178049693497

    return sum(p_r + p_g + p_b)
end

function powerR(R...)
    p_r = 0
    for r in R
        p_r += 9.33773772e-6*r^3 - 1.77241480e-3*r^2 + 3.63591812e-1*r + 0.40109990900690207
    end

    return p_r
end

function powerG(G...)
    p_g = 0
    for g in G
        p_g += 3.31277854e-6*g^3 + 9.59555435e-4*g^2 + 2.57421833e-2*g + 2.9828252112477784
    end

    return p_g
end

function powerB(B...)
    p_b = 0
    for b in B
        p_b += 4.19953259e-6*b^3 + 9.91716987e-4*b^2 + 1.99992989e-1*b + 7.100178049693497
    end

    return p_b
end


"""
Generates a reusable distance model from the provided point to any other point
within a given space.

Params:
    R: The red intensity of the color.
    G: The green instensity of the color.
    B: The blue intensity of the color.
    space: The space function that will represent the model used.
    metric: The metric used as a distance if non-custom models are used.
"""
function generate_distance(R::Float64, G::Float64, B::Float64, space, metric=nothing)
    # Check the metric used
    if metric == nothing
        # Obtain the point v coordinates in the space
        x, y, z = space(R, G, B)

        # Make the closure performant
        let x = x, y = y, z = z, space = space
            return function dist(R_n, G_n, B_n)
                        #=
                            Evaluates the euclidean distance between an RGB point
                            and the original point used to generate this function in
                            the required color space.

                            Params:
                                R_n: The red intensity of the color.
                                G_n: The green instensity of the color.
                                B_n: The blue intensity of the color.
                        =#

                        # Obtain the point position of the color in the required space
                        x_n, y_n, z_n = space(R_n, G_n, B_n)

                        # Calculate ΔE
                        return (x - x_n)^2 + (y - y_n)^2 + (z - z_n)^2
            end
        end
    else
        # Obtain the point v coordinates in the correct format
        original = RGB(R/255, G/255, B/255)

        # Make the closure performant
        let original = original, space = space
            return function dist(R_n, G_n, B_n)
                        #=
                            Evaluates the euclidean distance between an RGB point
                            and the original point used to generate this function in
                            the required color space.

                            Params:
                                R_n: The red intensity of the color.
                                G_n: The green instensity of the color.
                                B_n: The blue intensity of the color.
                        =#

                        # Obtain the point position of the color in the required space
                        new = RGB(R_n/255, G_n/255, B_n/255)

                        # Calculate ΔE
                        return Colors.colordiff(new, original; metric=space)
            end
        end
    end
end

"""
Generates a reusable parametrized distance model from the provided point to any
other point within a given space.

Params:
    R: The red intensity of the color.
    G: The green instensity of the color.
    B: The blue intensity of the color.
    space: The function that will convert the point to the given color space.
"""
function generate_parametrized_distance(R, G, B, space)
    # Obtain the point v coordinates in the space
    x, y, z = space(R, G, B)

    # Make the closure performant
    let x = x, y = y, z = z, space = space
        return function dist(R_n, G_n, B_n)
                    #=
                        Evaluates the euclidean distance between an RGB point
                        and the original point used to generate this function in
                        the required color space.

                        Params:
                            R_n: The red intensity of the color.
                            G_n: The green instensity of the color.
                            B_n: The blue intensity of the color.
                    =#

                    # Obtain the point position of the color in the required space
                    x_n, y_n, z_n = space(R_n, G_n, B_n)

                    # Calculate ΔE
                    return (x - x_n)^2 + (y - y_n)^2 + (z - z_n)^2
        end
    end
end


"""
A log barrier cost model.

Params:
    R: The red intensity of the color.
    G: The green instensity of the color.
    B: The blue intensity of the color.
    delta_E: The color difference as measured by some distance metric.
    eps: The maximum allowable deviation in delta_E.
    t: The factor that multiplies the power model and should tend to infinity.
"""
function log_cost(R, G, B, delta_E, eps, t=1e4)
    # If the constraint are violated we output infinity
    if eps < delta_E
        return Inf

    # Otherwise output normal power
    else
        return t * power(R, G, B) - log(eps - delta_E)
    end
end


"""
An exponential penalty cost model.

Params:
    R: The red intensity of the color.
    G: The green instensity of the color.
    B: The blue intensity of the color.
    delta_E: The color difference as measured by some distance metric.
    eps: The maximum allowable deviation in delta_E.
    t: The factor that multiplies the power model and should tend to infinity.
"""
function exp_cost(R, G, B, delta_E, eps, t=1e4)
    return t * power(R, G, B) + exp(max(0, delta_E-eps))
end


"""
A function that combines the costs and distance function and gives a value for
any pixel.
"""
# TODO finish the documentation
function generate_cost(model, distance, eps, metric)
    let model = model, distance = distance, eps = eps, metric = metric
        return function cost(RGB_pixel)
                    # If any color parts are negative return infinite cost as
                    # otherwise color conversions fail
                    if RGB_pixel[1] < 0 || RGB_pixel[2] < 0 || RGB_pixel[3] < 0
                        return Inf
                    end

                    # Otherwise conduct normal power calculation
                    delta_E = distance(RGB_pixel[1], RGB_pixel[2], RGB_pixel[3])
                    if metric == nothing
                        delta_E = 1.41 * delta_E^(0.63 * 0.5)
                    else
                        delta_E = delta_E
                    end
                    return model(RGB_pixel[1], RGB_pixel[2], RGB_pixel[3], delta_E, eps)
               end
    end
end
