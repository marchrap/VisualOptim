"""
An inline function that converts the RGB values to their linear counterparts

Params:
    V: the original RGB vector.
"""
@inline function rgb2linear(V)
    v = V / 255.

    if v <= 0.04045
        v *= 25.0 / 323.0
    else
        v = ((200.0 * v + 11.0) / 211.0)^2.4
    end

    return v
end


"""
Function that generates the function resposible for converting a point to
he CIECAM16 UCS space.
"""
function CAM16Model()
    # Define model parameters
    X_w = [95.047; 100.0; 108.883]
    Y_b = 20.0
    X_wr = Y_wr = Z_wr = 100.0
    L_A = 40.0
    F = 1.0
    c = 0.69
    N_c = 1.0
    M_16 = [0.401288 0.650173 -0.051461;
            -0.250268 1.204414 0.045854;
            -0.002079 0.048952 0.953127]

    # Calculate all values independent of the input
    R_w = M_16 * X_w
    D = F * (1 - 1 / 3.6 * exp((-L_A - 42) / 92))
    if D > 1
        D = 1
    elseif D < 0
        D = 0
    end
    D_R = D * X_w[2] ./ R_w .+ 1.0 .- D
    k = 1 / (5 * L_A + 1)
    F_L = k^4 * L_A + 0.1 * (1 - k^4)^2 * (5 * L_A)^(1 / 3)
    n = Y_b / X_w[2]
    z_w = 1.48 + sqrt(n)
    N_bb = 0.725 / (n^0.2)
    N_cb = N_bb
    R_wc = D_R .* R_w
    R_aw = 400.0 .* ((F_L .* R_wc) ./ 100.0).^0.42 ./ (((F_L .* R_wc) ./ 100.0).^0.42 .+ 27.13)
    A_w = (2 * R_aw[1] + R_aw[2] + 1 / 20 * R_aw[3]) * N_bb

    let D_G = D_R[2], D_B = D_R[3], D_R = D_R[1], n = n, F_L = F_L, z_w = z_w,
          c = c, N_bb = N_bb, A_w = A_w

        return  function rgb2cam(V_r, V_g, V_b)
                    #=
                        Converts the RGB values to CIECAM16 UCS values.

                        Params:
                            V_r: The red component of the color.
                            V_g: The green component of the color.
                            V_b: The blue component of the color.
                    =#
                    # Move RGB values to linear ones
                    r = rgb2linear(V_r)
                    g = rgb2linear(V_g)
                    b = rgb2linear(V_b)

                    # Convert linear RBG values to XYZ
                    x = 0.4124564 * r + 0.3575761 * g + 0.1804375 * b
                    y = 0.2126729 * r + 0.7151522 * g + 0.0721750 * b
                    z = 0.0193339 * r + 0.1191920 * g + 0.9503041 * b

                    # TODO change namign to be correct with the paper
                    # Calculate cone responses
                    R_y = 0.401288 * x + 0.650173 * y + -0.051461 * z
                    G_y = -0.250268 * x + 1.204414 * y + 0.045854 * z
                    B_y = -0.002079 * x + 0.048952 * y + 0.953127 * z

                    # Color adaptation
                    R_c = D_R * R_y
                    G_c = D_G * G_y
                    B_c = D_B * B_y

                    # Postadaptation
                    R_a = 400.0 * sign(R_c) * (F_L * abs(R_c))^0.42 / ((F_L * abs(R_c))^0.42 + 27.13)
                    G_a = 400.0 * sign(G_c) * (F_L * abs(G_c))^0.42 / ((F_L * abs(G_c))^0.42 + 27.13)
                    B_a = 400.0 * sign(B_c) * (F_L * abs(B_c))^0.42 / ((F_L * abs(B_c))^0.42 + 27.13)

                    # Components
                    a = R_a - 12.0 / 11.0 * G_a + 1.0 / 11.0 * B_a
                    b = 1.0 / 9.0 * R_a + 1.0 / 9.0 * G_a - 2.0 / 9.0 * B_a
                    h = atan(b, a)

                    # Eccentricity
                    e_t = 1.0 / 4.0 * (cos(h + 2.0) + 3.8)

                    # Chromatic resposne
                    A = (2.0 * R_a + G_a + 1.0 / 20.0 * B_a) * N_bb

                    # Lightness
                    J = 100.0 * (A / A_w)^(c * z_w)

                    # Colorfullness
                    t = 50000.0 / 13.0 * N_c * N_bb * e_t * sqrt(a^2 + b^2) / (R_a + G_a + 21.0 / 20.0 * B_a + 0.305)
                    α = t^0.9 * (1.64 - 0.29^n)^0.73
                    C1 = α * sqrt(J / 100.0)
                    M = C1 * F_L^0.25

                    # Move to UCS
                    J_prime = 1.7 * J / (1 + 0.007 * J)
                    M_prime = log(1 + 0.0228 * M) / 0.0228
                    a_prime = M_prime * cos(h)
                    b_prime = M_prime * sin(h)

                    # Return the result
                    return J_prime, a_prime, b_prime
        end
    end
end
