#include <gtest/gtest.h>
#include "../src/Compressed_spars_row/CompressedMatrix.h"
#include "../src/Compressed_spars_row/CompressedSorting.h"


TEST(CompressedTest, firstTest){
    std::vector<element<int>> matrix = {{0,0,1}, {0,1,2},
                                        {1,1,1}, {1,2,3},
                                        {2,0,4}, {2,2,2}};
    sort_me_plz(matrix);
    CompressedMatrix<int> res (matrix, 3, 3);

    std::vector<float> vector_of_values = {1, 2, 1, 3, 4, 2};
    std::vector<unsigned int> colomn_index_vector = {0, 1, 1, 2, 0, 2};
    std::vector<unsigned int> counts_vector = {0, 2, 4, 6};

    for (int i = 0; i < vector_of_values.size(); i++){
        ASSERT_NEAR(vector_of_values[i], res.getValue()[i], 0.001);
        ASSERT_NEAR(colomn_index_vector[i], res.get_index()[i], 0.001);
    }
    for(int i = 0; i < counts_vector.size(); ++i){
        ASSERT_NEAR(counts_vector[i], res.get_count()[i], 0.001);
    }
}

TEST(CompressedTest, secondTest){
    std::vector<element<float>> matrix = {{0,0,1}, {0,1,2},
                                        {1,1,1}, {1,2,3},
                                        {2,0,4}, {2,2,2}};
    sort_me_plz(matrix);
    CompressedMatrix<float> res (matrix, 3, 3);

    std::vector<float> vector_of_values = {1, 2, 1, 3, 4, 2};
    std::vector<unsigned int> colomn_index_vector = {0, 1, 1, 2, 0, 2};
    std::vector<unsigned int> counts_vector = {0, 2, 4, 6};
    std::vector<float> x = {1, 3, 5};
    std::vector<float> solve = {7, 8, 14};
    std::vector<float> result = res.dot(x);

    for (int i = 0; i < vector_of_values.size(); i++){
        ASSERT_NEAR(vector_of_values[i], res.getValue()[i], 0.001);
        ASSERT_NEAR(colomn_index_vector[i], res.get_index()[i], 0.001);
    }
    for(int i = 0; i < counts_vector.size(); ++i){
        ASSERT_NEAR(counts_vector[i], res.get_count()[i], 0.001);
    }

}


TEST(CompressedTest, thirdTest){
    std::vector<element<float>> matrix = {{0,0,15}, {0,1,1.2}, {0,4,1},
                                          {1,3,1}, {1,4,2},
                                          {2,1,3}, {2,2,4}, {2,4,7},
                                          {3,0,1}, {3,1,2}, {3,3,1},
                                          {4,0,5.8}, {4,1,3},{4,4,9}
    };
    sort_me_plz(matrix);
    CompressedMatrix res(matrix, 5, 5);


    std::vector<float> vector_of_values = {15, 1.2, 1, 1, 2, 3, 4, 7, 1, 2, 1, 5.8, 3, 9};
    std::vector<unsigned int> colomn_index_vector = {0, 1, 4, 3, 4, 1, 2, 4, 0, 1, 3, 0, 1, 4};
    std::vector<unsigned int> counts_vector = {0, 3, 5, 8, 11, 14};
    for (int i = 0; i < vector_of_values.size(); i++){
        ASSERT_NEAR(vector_of_values[i], res.getValue()[i], 0.001);
        ASSERT_NEAR(colomn_index_vector[i], res.get_index()[i], 0.001);
    }
    for(int i = 0; i < counts_vector.size(); ++i){
        ASSERT_NEAR(counts_vector[i], res.get_count()[i], 0.001);
    }
}
