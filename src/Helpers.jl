"""
Add helper functions such as:
- showing original image
- showing changed image
x plotting them next to each other
x creating histograms of their colors
- showing the change of lightness
- saving the images (maybe a struct with orignal and changed image)
x evaluating power gain per pixel
x evaluating power gain over all image
x evaluating power over the whole space
- saving an optimization to a file
- maybe a struct for the optimization space
x function optimizing the space
"""
# TODO Finish documentation
function optimize_space!(space, eps, color_model)
    # Create the progress bar
    progress_bar = Progress(length(space), 1, "Computing the matrix... ", 50)
    update!(progress_bar, 0)

    # Create the lock so that we can update the progress bar
    atomic = Threads.Atomic{Int}(0)
    atomic2 = Threads.Atomic{Int}(0)
    lock = Threads.SpinLock()

    # Run over the points
    # TODO split into multiple sub tasks rather than this as this has to read
    #      the actual memory!
    Threads.@threads for pixel in space
        # Define the distance and the cost
        distance = generate_distance(pixel.r * 255., pixel.g * 255., pixel. b * 255., color_model)
        cost = generate_cost(log_cost, distance, eps)

        # Obtain the solution and alter the space
        res = optimize(cost, [pixel.r * 255., pixel.g * 255., pixel.b * 255.], Optim.Options(g_tol = 1e-3))
        if !Optim.converged(res)
            Threads.atomic_add!(atomic2, 1)
        end
        solution = Optim.minimizer(res) ./ 255.
        space[convert(Int64, 255. * pixel.r) + 1, convert(Int64, 255 * pixel.g) + 1, convert(Int64, 255 * pixel.b) + 1] = RGB(solution[1], solution[2], solution[3])

        # Increase the atomic counter
        Threads.atomic_add!(atomic, 1)

        # Grab the lock and update the progress bar
        Threads.lock(lock)
        update!(progress_bar, atomic[])
        Threads.unlock(lock)
    end
    println("$(atomic2[]) out of $(length(space)) ($(atomic2[]/length(space))%) did not converge")
end

"""
A helper function that allows to read the the result from the optimized
array and protects against the NaN values.

Params:
    pixel: The RGB triplet that specifies which color is being optimized.
    optimized: The optimized array.
"""
function convert_pixel(pixel, optimized)
    # Obtain the value
    value = optimized[convert(Int64, 255*pixel.r)+1, convert(Int64, 255*pixel.g)+1, convert(Int64, 255*pixel.b)+1]

    # Check if the value is NaN
    if isnan(value.r)
        println("Found a NaN...")
        return pixel
    end
    return value
end


"""
Converts the image to its optimized version provided the optimized space and the
original image.
"""
function convert_image(image, optimized)
    return map(pixel -> convert_pixel(pixel, optimized), image)
end


"""
Compares the power of two images.

Params:
    img1: The first image as loaded by Julia.
    img2: The second image as loaded by Julia.
"""
function compare_power(img1, img2)
    original_power = 0
    new_power = 0

    # Evaluate the power
    for i in 1:size(img1)[1], j in 1:size(img1)[2]
        original_power += power(img1[i, j].r*255., img1[i, j].g*255., img1[i, j].b*255.)
        new_power += power(img2[i, j].r*255., img2[i, j].g*255., img2[i, j].b*255.)
    end

    # Print the results
    println("Original power $original_power; New power $new_power")

    return (new_power-original_power)/original_power
end


"""
Evaluates the power gains of the whole space instead of considering single
points.

Params:
    optimized: The optimized array.
"""
function evaluate_space_power(optimized)
    # Define the original and modified powers
    original = 0
    modified = 0

    # Loop over the spaces and evvaluate the power
    for i in 0:255, j in 0:255, k in 0:255
        original += power(i, j, k)
        pixel = optimized[i+1, j+1, k+1]
        modified += power(pixel.r * 255., pixel.g * 255., pixel.b * 255.)
    end

    return original, modified
end


"""
Evaluates the changes in the RGB compononents of a whole space.

Params:
    optimized: The optimized array.
"""
function evaluate_space_RGB_change(optimized)
    # Define the changes array and set the corresponding values to the required
    # RGB changes
    changes = zeros(256 * 256 * 256, 3)
    for i = 0:255, j = 0:255, k = 0:255
        pixel = optimized[i+1, j+1, k+1]
        changes[i*256^2 + j*256 + k+1, 1] = pixel.r * 255. - i
        changes[i*256^2 + j*256 + k+1, 2] = pixel.g * 255. - j
        changes[i*256^2 + j*256 + k+1, 3] = pixel.b * 255. - k
    end

    return changes
end


"""
Evaluates the histogram of an image and its colors and evaluates its other
characteristics.
"""
# TODO include things such as the lightness spectrum etc
function histogram(image)
    # Define the three arrays for the colors
    R = []
    G = []
    B = []

    # Iterate over all pixels
    for pixel in image
        append!(R, pixel.r * 255.)
        append!(G, pixel.g * 255.)
        append!(B, pixel.b * 255.)
    end

    return R, G, B
end


"""
Evaluates the relative change in RGB color values.
"""
function compare_colors(original_image, optimized_image)
    # Obtain the RGB histograms for both the images
    R, G, B = histogram(original_image)
    R2, G2, B2 = histogram(optimized_image)

    # Calculate and return the changes of the means
    delta_R = (mean(R) - mean(R2))/mean(R)
    delta_G = (mean(G) - mean(G2))/mean(G)
    delta_B = (mean(B) - mean(B2))/mean(B)

    return delta_R, delta_G, delta_B
end


"""
Plots the image together with its corresponding histograms.
"""
function image_histogram(image)
    # Obtain the R, G, B component lists
    R, G, B = histogram(img)

    # Plot the image and the histograms on the right hand side
    display(plot(plot(img),
        Plots.histogram(R, color="red", label="R"),
        Plots.histogram(G, color="green", label="G"),
        Plots.histogram(B, color="blue", label="B"),
        layout = @layout([ a [b; c; d]])))
end


"""
Plots the two images side by side.
"""
function plot_side_by_side(image1, image2)
    display(plot(plot(image1), plot(image2), layout = @layout([ a b ])))
end


"""
Compares two images by plotting their difference and histograms of their changes.
"""
# TODO check what happens if you plot the difference to be the CIECAM16 difference
# instead of the simple euclidean
function image_histogram(original_image, optimized_image)
    # Obtain the R, G, B component lists
    R, G, B = histogram(original_image)
    R_opt, G_opt, B_opt = histogram(optimized_image)

    # Obtain the difference image
    difference = zeros(size(original_image))
    for i in eachindex(difference)
        a = [original_image[i].r, original_image[i].g, original_image[i].b]
        b = [optimized_image[i].r, optimized_image[i].g, optimized_image[i].b]
        difference[i] = norm(a-b)*255
    end

    # Plot the difference image and the histograms on the right hand side
    display(plot(Plots.heatmap(difference),
        Plots.histogram(Any[R, R_opt], alpha=0.2, title="R"),
        Plots.histogram(Any[G, G_opt], alpha=0.2, title="G"),
        Plots.histogram(Any[B, B_opt], alpah=0.2, title="B"),
        layout = @layout([ a [b; c; d]])))
end


"""
Compares two images by plotting their difference and histograms of their changes
but only in a single color.
"""
# TODO check what happens if you plot the difference to be the CIECAM16 difference
# instead of the simple euclidean
function image_histogram_color(original_image, optimized_image, color)
    # Obtain the R, G, B component lists
    R, G, B = histogram(original_image)
    R_opt, G_opt, B_opt = histogram(optimized_image)

    # Obtain the difference image
    difference = zeros(size(original_image))
    for i in eachindex(difference)
        a = [original_image[i].r, original_image[i].g, original_image[i].b]
        b = [optimized_image[i].r, optimized_image[i].g, optimized_image[i].b]
        if color == "R":
            difference[i] = norm(a[1]-b[1])*255
        elseif color == "G"
            difference[i] = norm(a[2]-b[2])*255
        elseif color == "B"
            difference[i] = norm(a[2]-b[2])*255
        end
    end

    # Plot the difference image and the histograms on the right hand side
    display(plot(Plots.heatmap(difference),
        Plots.histogram(Any[R, R_opt], alpha=0.2, title="R"),
        Plots.histogram(Any[G, G_opt], alpha=0.2, title="G"),
        Plots.histogram(Any[B, B_opt], alpah=0.2, title="B"),
        layout = @layout([ a [b; c; d]])))
end
