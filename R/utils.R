
var_components <- function(values = NULL, q = 0.2, 
                           a = NULL, d = NULL) {
    p <- 1-q
    
    # Assigned values
    if(is.null(a) | is.null(d)) {
        midpoint <- (values[1] + values[3]) / 2
        a <- values[1] - midpoint
        d <- values[2] - midpoint
    }
    
    # VA
    alpha <- as.numeric(a + d * (q-p))
    va <- 2 * p * q * alpha^2
    
    # VD
    vd <- (2 * p * q * d)^2
    
    # VG
    vg <- va + vd
    result_df <- data.frame(
        Componentes = c("VG", "VA", "VD"),
        Valores = c(vg, va, vd)
    )
    return(result_df)
}

avg_effect <- function(a = NULL, d = NULL, q = NULL) {
    p <- 1 - q
    alpha <- as.numeric(a + d * (q-p))
    alpha_1 <- as.numeric(q * alpha)
    alpha_2 <- as.numeric(-p * alpha)
    df <- data.frame(alpha = alpha, 
                     alpha1 = alpha_1, 
                     alpha2 = alpha_2)
    return(df)
}

# G, A e D
calc_GAD <- function(a = NULL, d = NULL, q = NULL, alpha = NULL) {
    
    p <- 1 - q
    # Breeding values for each genotype
    A <- c(a1a1 = 2*q*alpha, 
           a1a2 = alpha * (q-p),
           a2a2 = -2 * p * alpha)
    
    # Genotypic value in terms of deviation
    G <- c(a1a1 = 2*q*(alpha - q*d),
           a1a2 = alpha*(q-p) + 2*p*q*d,
           a2a2 = -2*p*(alpha + p*d))
    
    # Dominance deviation
    D <- c(a1a1 = -2 * q^2 * d,
           a1a2 = 2*p*q*d,
           a2a2 = -2 * p^2 * d)
    
    df <- data.frame(
        G = G,
        A = A,
        D = D
    )
    return(df)
}


means_and_values <- function(genotype_values = NULL, p = NULL) {
    if(length(genotype_values) != 3) {
        stop("Object 'genotype_values' must be a character vector of length 3.")
    }
    q <- 1 - p
    
    # Assigned values
    midpoint <- (genotype_values[1] + genotype_values[3]) / 2
    a <- genotype_values[1] - midpoint
    d <- genotype_values[2] - midpoint
    
    # Mean
    M <- a * (p-q) + 2 * d * p * q
    M_plus_mp <- as.numeric(midpoint + M)
    Mean <- c(Mean = as.numeric(M),
              Mean_plus_midpoint = M_plus_mp)
    
    # Degree of dominance  
    dg <- as.numeric(d / a)
    
    # Average effect
    alpha <- as.numeric(a + d * (q-p))
    alpha1 <- as.numeric(q * alpha)
    alpha2 <- as.numeric(-p * alpha)
    
    # Breeding values for each genotype
    A <- c(a1a1 = 2*q*alpha, 
           a1a2 = alpha * (q-p),
           a2a2 = -2 * p * alpha)
    
    # Genotypic value in terms of deviation
    G <- c(a1a1 = 2*q*(alpha - q*d),
           a1a2 = alpha*(q-p) + 2*p*q*d,
           a2a2 = -2*p*(alpha + p*d))
    
    # Dominance deviation
    D <- c(a1a1 = -2 * q^2 * d,
           a1a2 = 2*p*q*d,
           a2a2 = -2 * p^2 * d)
    
    table <- data.frame(
        Genotype = c("A1A1", "A1A2", "A2A2"),
        Assigned_values = c(paste0("+a = ", a),
                            paste0("d = ", d),
                            paste0("-a = ", -a)),
        A = A,
        G = G,
        D = D
    )
    rownames(table) <- NULL
    result_list <- list(
        Mean = Mean,
        `Degree of dominance` = dg,
        `Average effects` = c(alpha = alpha,
                              alpha1 = alpha1,
                              alpha2 = alpha2),
        Summary_table = table
    )
    return(result_list)
}

