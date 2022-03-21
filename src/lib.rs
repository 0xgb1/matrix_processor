fn transpose(M: &Vec<Vec<i64>>) -> Vec<Vec<i64>> {
    let mut matrix = def_2D_matrix(M[0].len(), M.len());
    let mut count = 0;
    for i in 0..M.len() {
        for j in 0..M[0].len() {
            if i == j {
                matrix[i][i] = M[i][i];
                continue;
            }
            else {
                matrix[j][i] = M[i][j];
            }
        }
        count += 1;
    }
    matrix
}

fn dot_product(v1: &Vec<i64>, v2: &Vec<i64>) -> i64 {
    if v1.len() != v2.len(){
        panic!("These vectors are incompatible for dot product");
    }
    //Error to be made better later, dependent on GUI
   
    else { 
        let mut result: i64 = 0;
        for i in 0..v1.len() {
            result += v1[i] * v2[i];
        }
        result
    }
}

pub fn def_2D_matrix(m: usize, n: usize) -> Vec<Vec<i64>> {

    //Make Error Case for either m or n being equal to 0

    let mut A = Vec::with_capacity(m);
    for _index in 0..m {
        A.push(vec![0_i64; n]);
    }
    A
}

pub fn multiply_matrices(A: Vec<Vec<i64>>, B: &Vec<Vec<i64>>) -> Vec<Vec<i64>> {
    if A[0].len() != B.len() {
        panic!("These matrices are not able to be multiplied!");
    }
    else {
        //let mut new_B = def_2D_matrix(B.len(), B[0].len());
        let b_ref = &*B;
        let new_B = transpose(b_ref);
        let mut result = def_2D_matrix(A.len(), new_B.len());
        for i in 0..A.len() {
            for j in 0..new_B.len() {
                result[i][j] = dot_product(&A[i], &new_B[j]);
            }
        }
        result
    }
}

pub fn matrix_power(A: Vec<Vec<i64>>, n: usize) -> Vec<Vec<i64>> {
    if A.len() != A[0].len(){
        panic!("This matrix cannot be raised to a power!");
    }
    else {
        let mut result = A.clone();
        for _i in 0..n-1 {
            let tmp = &A;
            result = multiply_matrices(result, tmp);
        }
        result
    }
}

pub fn determinant(A: Vec<Vec<i64>>) -> i64 {
    if A.len() != A[0].len() {
        panic!("This matrix isn't square!");
    }
    let mut result: i64 = 0;
    if A.len() == 2 {
        result = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    }
    else if A.len() > 2 {
        let mut res: i64 = 0;
        let mut multiplier = 1;
        for i in 0..A.len() {
            let first = A[0][i];
            let mut mini_det_vec = vec![];
                for j in 1..A.len(){
                    let mut first_vec = vec![];
                    for k in 0..i {
                        first_vec.push(A[j][k]);
                    }
                    if i < A.len() - 1 {
                        for k in i+1..A.len(){
                            first_vec.push(A[j][k]);
                        }
                    }
                    mini_det_vec.push(first_vec);
                }
            res += first * determinant(mini_det_vec) * multiplier;
            multiplier *= -1;
            result = res;
        }
    }
    result
}

pub fn matrix_inverse(A: Vec<Vec<i64>>) -> (Vec<Vec<i64>>, u8, i64) {
    if A.len() != A[0].len() || A.len() == 1{
        panic!("This matrix isn't square or it is too small!");
    }
    else {
        let mut new_matrix = A.clone();
        for i in 0..A.len() {
            for j in 0..A.len(){
                let mut minor = vec![];
                for k in 0..i {
                    let mut minor_row = vec![];
                    for l in 0..j {
                        minor_row.push(A[k][l]);
                    }
                    if j < A.len() - 1 {
                        for l in j+1..A.len(){
                            minor_row.push(A[k][l]);
                        }
                    }
                    minor.push(minor_row);
                }
                if i < A.len() - 1 {
                    for k in i+1..A.len() {
                        let mut minor_row = vec![];
                        for l in 0..j {
                            minor_row.push(A[k][l]);
                        }
                        if j < A.len() - 1 {
                            for l in j+1..A.len(){
                                minor_row.push(A[k][l]);
                            }
                        }
                        minor.push(minor_row);
                    }
                }
            new_matrix[i][j] = determinant(minor);
            }
        }
        let mut multiplier: i8 = 1;
        for i in 0..A.len() {
            for j in 0..A.len(){
                new_matrix[i][j] *= multiplier as i64;
                multiplier *= -1;
            }
            if A.len() % 2 == 0 {
                multiplier *= -1;
            }
        }
        new_matrix = transpose(&new_matrix);
        (new_matrix, 1, determinant(A))
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn multiply_matrices_test() {
        assert_eq!(super::multiply_matrices(vec![vec![5, 6],vec![1, 8],vec![3, 7],vec![4, 2]], 
            &vec![vec![5, 1, 3, 4], vec![6, 8, 7, 2]]),
            vec![vec![61, 53, 57, 32], vec![53, 65, 59, 20], vec![57, 59, 58, 26], vec![32, 20, 26, 20]])
    }

    #[test]
    fn matrix_power_test() {
        assert_eq!(super::matrix_power(vec![vec![1, 0, 5, 0], vec![0, 4, 4, 9], 
            vec![0, 6, 9, 0], vec![1, 0, 4, 0]], 12),
            vec![vec![291867098311, 5409700623330, 8497849868540, 3626589826980], 
            vec![534933238905, 9914759068018, 15574687219726, 6646645589472], 
            vec![725317965396, 13443240216246, 21117459248223, 9012143156598], 
            vec![237841556797, 4408351383708, 6924875992247, 2955297403098]])
    }
    #[test]
    fn determinant_test() {
        assert_eq!(super::determinant(vec![vec![1,4,5,0,3,7], vec![4,0,5,2,0,4], 
            vec![3,9,1,0,0,7], vec![3,1,5,0,4,9], 
            vec![8,0,8,2,7,10], vec![1, 3, 4, 6, 0, 8]]),
            -25928)
    }
    #[test]
    fn matrix_inverse_test(){
        assert_eq!(super::matrix_inverse(vec![vec![4, 8, 1, 3], vec![5, 1, 0, 6], 
            vec![7, 3, 2, 4], vec![2, 0, 9, 1]]),
            (vec![vec![128, 194, -406, 76], vec![-184, 42, 74, 4], 
            vec![-20, -2, 54, -144], vec![-76, -370, 326, -64]], 1, -1208)
        )
    }
}