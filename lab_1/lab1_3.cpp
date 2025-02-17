#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>

using namespace std;

// Функция для генерации массива случайных чисел в интервале от -1 до 1
vector<double> generate_random_array(int size) {
    random_device rd;
    mt19937 engine(time(0));
    uniform_real_distribution<double> gen(-1.0, 1.0);

    vector<double> array(size);
    for (int i = 0; i < size; ++i) {
        array[i] = gen(engine);
    }
    return array;
}

// Функция сортировки вставками с измерением количества обменов и проходов
void insertion_sort(vector<double>& array, size_t& swaps, size_t& passes) {
    swaps = 0;
    passes = 0;
    for (size_t i = 1; i < array.size(); i++) {
        double key = array[i];
        int j = i - 1;
        
      if (array[j] > key) {
            swaps++;
        }  

        while (j >= 0 && array[j] > key) {
            array[j + 1] = array[j];
            j--;
      //      swaps++;
        }
        array[j + 1] = key;
        passes++;
    }
}

int main() {
    int attempts_per_series = 1;
    vector<int> array_sizes = { 1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000 };

    // Сохранение данных в файл для построения графиков
    ofstream time_file("times.txt");
    ofstream swap_file("swaps.txt");
    ofstream pass_file("passes.txt");

    for (int i = 0; i < size(array_sizes); i++) {
        size_t array_size = array_sizes[i];

        for (int attempt = 0; attempt < attempts_per_series; attempt++) {
            auto array = generate_random_array(array_size);

            cout << array_size << endl;

            chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
            size_t swap_count = 0, pass_count = 0;
            insertion_sort(array, swap_count, pass_count);
            chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();

            chrono::duration<double> sec_diff = end - start;

            time_file << array_sizes[i] << " " << sec_diff.count() << "\n";
            swap_file << array_sizes[i] << " " << swap_count << "\n";
            pass_file << array_sizes[i] << " " << pass_count << "\n";

            cout << sec_diff.count() << endl;
        }
    }

    time_file.close();
    swap_file.close();
    pass_file.close();
    return 0;
};
