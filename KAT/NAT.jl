using Symbolics
using LinearAlgebra
using SHA
using Statistics

message = "Hello"
@variables x[1:64]

function modify_coeffcients(A, q)
    coeffs_map = Symbolics.value(A).dict#turn the polynomial into matrix
    B = 0#create our new polynomial
    for (key, value) in coeffs_map#iterate through the dict(monomial, coeffcients)
        new_coef = value % q
        B += new_coef * key
    end
    return B
end

function sparse_polynomial(n, t, b, q)
    @variables x[1:n]
    sparse_poly = 0

    for i in 1:t #we want t monomials
        deg = rand(1:b)

        monomial_d = 1
        for j in 1:deg
            monomial_d *= x[rand(1:n)]
        end
        sparse_poly += rand(1:(q-1)) * monomial_d

    end
    #     println(sparse_poly)
    A = modify_coeffcients(sparse_poly, q)
    return A
end

function num_identity(k)
    a = zeros(Num, k)
    a[1] = 1
    for i in 2:k
        b = zeros(Num, k)
        b[i] = 1
        a = [a b]
    end
    return a
end


function generate_U(k, n, t, B, q)
    Es = []
    for i in 1:(k-1)
        for j in (i+1):k
            I_matrix = num_identity(k)
            I_matrix[i, j] = sparse_polynomial(n, t, B, q)
            push!(Es, I_matrix)
        end
    end
    return Es
end

function generate_L(k, n, t, B, q)
    Es = []
    for i in 2:k
        for j in 1:(i-1)
            I_matrix = num_identity(k)
            I_matrix[i, j] = sparse_polynomial(n, t, B, q)
            push!(Es, I_matrix)
        end
    end
    return Es
end


function cols_to_keep(k, l)
    col_opts = [i for i in 1:k]
    for i in 1:(k-l)
        random_number = size(col_opts, 1)
        splice!(col_opts, rand(1:random_number))
    end
    return col_opts
end

function random_del_cols(k, l, matrix)
    #     col_opts=[i for i in 1:l]
    col_opts = cols_to_keep(k, l)
    matrix_kl = zeros(Num, k, l)

    j = 1
    for i in col_opts
        matrix_kl[:, j] = matrix[:, i]
        j += 1
    end
    return col_opts, matrix_kl
end

function find_kl_inverse(to_keep, S_inv)
    k = size(S_inv, 1)#find the number of rows
    ml_inv = zeros(Num, length(to_keep), k)

    j = 1
    for i in to_keep
        ml_inv[j, :] = S_inv[i, :]
        j += 1
    end
    return ml_inv
end

function modify_neg_coeff(A, q)
    coeffs_map = Symbolics.value(A).dict#turn the polynomial into matrix
    B = 0#create our new polynomial
    for (key, value) in coeffs_map#iterate through the dict(monomial, coeffcients)
        if value > q || value < -q
            new_coef = mod(value, q)
            B += new_coef * key
        else
            B += key * value
        end
    end
    B
end

function modify_matrix_coef(matrix, q)
    res = simplify.(expand.(matrix))#expand and then simplify
    for i in 1:Int(size(res, 1))
        for j in 1:Int(size(res, 2))
            a = res[i, j]

            if isone(a) == false && iszero(a) == false
                res[i, j] = modify_neg_coeff(a, q)
            end
        end
    end
    return res
end



function matrix_kl(k, l, n, t, B, q)
    #generate upper unitriangular matrix
    U = generate_U(k, n, t, B, q)
    res_U = U[1]
    for i in 2:Int(size(U)[1])
        res_U = res_U * U[i]
    end

    #generate lower unitriangular matrix
    L = generate_L(k, n, t, B, q)
    res_L = L[1]
    for i in 2:Int(size(L)[1])
        res_L = res_L * L[i]
    end

    X = res_U * res_L


    uSize = Int(size(U)[1])
    invU = inv(U[uSize])
    for i in 1:(uSize-1)
        invU = invU * (inv(U[uSize-i]))
    end


    lSize = Int(size(L)[1])
    invL = inv(L[lSize])
    for i in 1:(lSize-1)
        invL = invL * (inv(L[lSize-i]))
    end


    #     println(simplify.(expand.(invU*res_U)))
    #     println(simplify.(expand.(invL*res_L)))
    Y = (invL) * (invU)
    #     println(simplify.(expand.(X*Y)))


    col_keep, kl_matrix = random_del_cols(k, l, X)
    kl_inverse = find_kl_inverse(col_keep, Y)
    return X, Y, kl_matrix, kl_inverse
end



function hashing512(message)
    hash = sha2_512(message)
    bitstring512 = ""
    for i in hash
        bitstring512 = bitstring512 * bitstring(i)
    end

    A512 = []
    for i in bitstring512
        append!(A512, Int(i) - 48)
    end
    return A512
end

function polys(message)
    coeffs = []
    A512 = hashing512(message)
    last12 = A512[501:512,]
    for i in 0:3
        b3 = last12[(1+3*i):3*(i+1),]
        coef = 0
        for j in 1:3
            coef = coef + 2^(3 - j) * b3[j]
        end
        coef = mod(coef, 6)
        append!(coeffs, coef)
    end
    #     println(coeffs)

    first300 = A512[1:300]
    vars = []
    for i in 1:50
        v = first300[(1+6*(i-1)):(6*i),]
        vn = 0
        for j in 1:6
            vn = vn + 2^(6 - j) * v[j]
        end
        vn = vn + 1
        append!(vars, vn)
    end
    #     println(vars)

    middle200 = A512[301:500]
    mons = []
    for i in 1:5
        b40 = middle200[(1+40*(i-1)):(40*i)]

        for j in 1:4
            b10 = b40[(1+10*(j-1)):(10*j)]
            mon = 1

            for k in 1:10
                if b10[k] == 1

                    if vars[10*(j-1)+k] != 0
                        mon = mon * x[vars[10*(j-1)+k]]
                    end
                end
            end
            append!(mons, mon)

        end
    end
    #     println(mons)

    polys = []
    for i in 1:5
        pol = 0
        for j in 1:4
            pol = pol + mons[4*(i-1)+j] * coeffs[j]
        end
        append!(polys, pol)
    end


    Un = zeros(Num, 1, 3)
    for i in 1:3
        Un[1, i] = polys[i]
    end

    #     Un=simplify.(expand.(Un))
    U = zeros(Num, 3, 1)
    for i in 1:3
        U[i, 1] = polys[i]
    end
    return U, Un
end


A, B, C, D = matrix_kl(5, 3, 64, 2, 3, 6)
U, Un = polys(message)
V = Un * D

file = open("/Users/jchen056/scrap/KAT/publicKey.txt", "w")
for i in 1:Int(size(C, 1))
    for j in 1:Int(size(C, 2))
        write(file, string(C[i, j]))
        write(file, "\n")
    end
    #     write(file,"\n")
end
close(file)
filesize("/Users/jchen056/scrap/KAT/publicKey.txt.txt")

file = open("/Users/jchen056/scrap/KAT/privateKey.txt", "w")
for i in 1:Int(size(D, 1))
    for j in 1:Int(size(D, 2))
        write(file, string(D[i, j]))
        write(file, "\n")
    end
    #     write(file,"\n")
end
close(file)
filesize("/Users/jchen056/scrap/KAT/privateKey.txt.txt")

function verification(C, V, Un)
    cnt = 0
    for i in 1:20
        substitution_dict = Dict()
        for i in 1:64
            substitution_dict[x[i]] = rand(0:1)
        end

        W = V * C
        left = substitute.(W, (substitution_dict,))
        right = substitute.(Un, (substitution_dict,))
        if isequal(left, right)
            cnt = cnt + 1
        end
    end
    return println(cnt == 20)
end

@time verification(C, V, Un)