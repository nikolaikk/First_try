//#include<conio.h>
#include<iostream>
#include <fstream>      // std::ifstream

using namespace std;

int main()

{
    ifstream file;
    file.open( "/home/rzlin/yk83ened/Project Nik/cpp/const_to_simplify_RE.txt");
    char test[5] = {0};

    for( int i=0; i<50; i++)
        file >> test[i];

    file.close();
    for( int i=0; i<50; i++)
        cout << test[i] << endl;


//    int x, y;
//    ifstream in("const_to_simplify_RE.txt");

//    if (!in) {
//    cout << "Cannot open file.\n";
//    return;
//    }

//    for (y = 0; y < 15; y++) {
//    for (x = 0; x < 15; x++) {
//      in >> distances[x][y];
//    }
//    }

//    in.close();

//    array_2d = new int*[200];

//    ifstream file('const_to_simplify_RE.txt');

//    for (unsigned int i = 0; i < 200; i++) {
//        array_2d[i] = new int[200];

//        for (unsigned int j = 0; j < 200; j++) {
//            file >> array_2d[i][j];
//        }
//    }


//    int a[10][10], b[10][10], c[10][10];

//    int x, y, i, j, m, n;

//    cout << "\nEnter the number of rows and columns for Matrix A:::\n\n";
//    cin >> x >> y;

//    // x denotes number rows in matrix A
//    // y denotes number columns in matrix A
//    cout << "\n\nEnter elements for Matrix A :::\n\n";
//    for (i = 0; i < x; i++)
//    {
//        for (j = 0; j < y; j++)
//        {
//            cin >> a[i][j];
//        }
//        cout << "\n";
//    }
//    cout << "\n\nMatrix A :\n\n";

//    for (i = 0; i < x; i++)
//    {
//        for (j = 0; j < y; j++)
//        {
//            cout << "\t" << a[i][j];
//        }
//        cout << "\n\n";
//    }

//    cout << "\n-----------------------------------------------------------\n";
//    cout << "\nEnter the number of rows and columns for Matrix B:::\n\n";
//    cin >> m >> n;

//    // m denotes number rows in matrix B
//    // n denotes number columns in matrix B
//    cout << "\n\nEnter elements for Matrix B :::\n\n";

//    for (i = 0; i < m; i++)
//    {
//        for (j = 0; j < n; j++)
//        {
//            cin >> b[i][j];
//        }
//        cout << "\n";
//    }

//    cout << "\n\nMatrix B :\n\n";
//    for (i = 0; i < m; i++)
//    {
//        for (j = 0; j < n; j++)
//        {
//            cout << "\t" << b[i][j];
//        }
//        cout << "\n\n";
//    }

//    if (y == m)
//    {
//        for (i = 0; i < x; i++)
//        {
//            for (j = 0; j < n; j++)
//            {
//                c[i][j] = 0;
//                for (int k = 0; k < m; k++)
//                {
//                    c[i][j] = c[i][j] + a[i][k] * b[k][j];
//                }
//            }
//        }

//        cout << "\n-----------------------------------------------------------\n";
//        cout << "\n\nMultiplication of Matrix A and Matrix B :\n\n";

//        for (i = 0; i < x; i++)
//        {
//            for (j = 0; j < n; j++)
//            {
//                cout << "\t" << c[i][j];
//            }
//            cout << "\n\n";
//        }
//    }

//    else
//    {
//        cout << "\n\nMultiplication is not possible";
//    }

    return 0;
}
