#include"Header.h"
#include<iostream>
#include<vector>
#include<algorithm>
#define MAX 300000

using namespace std;

//Selection Sort
void SelectionSort(int a[], int n) {
	for (int i = 0; i < n - 1; ++i) {
		int MIN = i;
		for (int j = i + 1; j < n; ++j)
			if (a[MIN] > a[j])
				MIN = j;
		swap(a[i], a[MIN]);
	}
}

//Insertion Sort
void InsertionSort(int a[], int n) {
	for (int i = 1; i < n; ++i) {
		int v = a[i], j = i - 1;
		while (j >= 0 && a[j] > v) {
			a[j + 1] = a[j];
			--j;
		}
		a[j + 1] = v;
	}
}

//Binary Insertion Sort
int BinarySearch(int a[], int n, int x) {
	int l = 0, r = n;
	while (l < r) {
		int m = (l + r) / 2;
		if (a[m] > x) {
			r = m;
			if (m == 0 || a[m - 1] <= x)
				return m;
		}
		else if (a[m] < x) {
			l = m + 1;
			if (m == n - 1 || a[m + 1] >= x)
				return m + 1;
		}
		else {
			return m;
		}
	}
	return n;
}
void BinaryInsertionSort(int a[], int n) {
	for (int i = 1; i < n; ++i) {
		int v = a[i], j = i - 1, end = BinarySearch(a, i, a[i]);
		for (; j >= end; --j)
			a[j + 1] = a[j];
		a[j + 1] = v;
	}
}

//Bubble Sort
void BubbleSort(int a[], int n) {
	for (int i = 1; i < n; ++i)
		for (int j = n - 1; j >= i; --j)
			if (a[j] < a[j - 1])
				swap(a[j], a[j - 1]);
}

//Shaker Sort
void ShakerSort(int a[], int n) {
	int i, j, k;
	for (i = 0; i < n;) {
		for (j = i + 1; j < n; ++j) {
			if (a[j] < a[j - 1])
				swap(a[j], a[j - 1]);
		}
		--n;
		for (k = n - 1; k > i; --k) {
			if (a[k] < a[k - 1])
				swap(a[k], a[k - 1]);
		}
		++i;
	}
}

//Shell Sort
void ShellSort(int a[], int n) {
	int h, i, j, temp;
	for (h = n / 2; h > 0; h /= 2) {
		for (i = h; i < n; ++i) {
			temp = a[i];
			for (j = i; j >= h && a[j - h] > temp; j -= h) {
				a[j] = a[j - h];
			}
			a[j] = temp;
		}
	}
}

//Heap Sort
void Sift(int a[], int l, int r) {
	int i = l, j = 2 * l + 1, x = a[l];
	while (j <= r) {
		if (j < r && a[j] < a[j + 1]) ++j;
		if (a[j] <= x) break;
		a[i] = a[j];
		i = j; j = i * 2;
	}
	a[i] = x;
}

void HeapSort(int a[], int n) {
	int l = n / 2 - 1;
	while (l >= 0) {
		Sift(a, l, n - 1);
		l--;
	}
	int r = n - 1;
	while (r > 0) {
		swap(a[0], a[r]);
		r--;
		Sift(a, 0, r);
	}
}

//Merge Sort
void MergeSort(int a[], int n) {
	if (n == 1 || !n) return;
	MergeSort(a, n / 2);
	MergeSort(a + n / 2, n - n / 2);
	int *b = new int[n], i = 0, j = n / 2, k = 0;
	while (i < n / 2 && j < n) {
		if (a[i] < a[j])
			b[k++] = a[i++];
		else
			b[k++] = a[j++];
	}
	while (i < n / 2)
		b[k++] = a[i++];
	while (j < n)
		b[k++] = a[j++];
	for (int i = 0; i < n; i++)
		a[i] = b[i];
	delete[]b;
}

//Quick Sort
void QuickSort(int a[], int l, int r) {
	if (l < r) {
		int i = l, j = r, pivot = a[(l + r) / 2];
		do {
			while (a[i] < pivot) ++i;
			while (a[j] > pivot) --j;
			if (i <= j)
				swap(a[i++], a[j--]);
		} while (i <= j);
		QuickSort(a, l, j);
		QuickSort(a, i, r);
	}
}

//Counting Sort
void CountingSort(int a[], int n) {
	int m = a[0];
	for (int i = 0; i < n; ++i)
		m = max(m, a[i]);
	vector<int>f(m + 1, 0), b(n);
	for (int i = 0; i < n; ++i)
		++f[a[i]];
	for (int i = 1; i <= m; ++i)
		f[i] += f[i - 1];
	for (int i = n - 1; i >= 0; --i) {
		b[f[a[i]] - 1] = a[i];
		--f[a[i]];
	}
	for (int i = 0; i < n; ++i)
		a[i] = b[i];
}

//Radix Sort
int GetDigit(int n, int d) {
	return (n / int(pow(10, d - 1))) % 10;
}

void LSD(int a[], int n, int d) {
	int f[10] = { 0 };
	vector<int>b(n);
	for (int i = 0; i < n; ++i)
		++f[GetDigit(a[i], d)];
	for (int i = 1; i <= 9; ++i)
		f[i] += f[i - 1];
	for (int i = n - 1; i >= 0; --i) {
		b[f[GetDigit(a[i], d)] - 1] = a[i];
		--f[GetDigit(a[i], d)];
	}
	for (int i = 0; i < n; ++i)
		a[i] = b[i];
}

void RadixSort(int a[], int n) {
	int m = a[0];
	for (int i = 1; i < n; ++i)
		m = max(m, a[i]);
	int d = log10(m) + 1;
	for (int i = 1; i <= d; ++i)
		LSD(a, n, i);
}

//Flash Sort
void FlashSort(int a[], int n) {
	int m = 0.43 * n;
	int* L = new int[m] {};
	int minA = a[0], maxA = a[0];
	for (int i = 1; i < n; ++i) {
		if (minA > a[i]) minA = a[i];
		if (maxA < a[i]) maxA = a[i];
	}
	
	for (int i = 0; i < n; ++i) {
		int k = 1ll * (m - 1) * (a[i] - minA) / (maxA - minA);
		++L[k];
	}
	for (int i = 1; i < m; ++i)
		L[i] += L[i - 1];
	
	int cnt = 0, i = 0, k = m - 1;
	while (cnt < n) {
		while (i > L[k] - 1) {
			++i;
			k = 1ll * (m - 1) * (a[i] - minA) / (maxA - minA);
		}
		int x = a[i]; //bat dau chu trinh
		while (i < L[k]) {
			k = 1ll * (m - 1) * (x - minA) / (maxA - minA);
			int y = a[L[k] - 1];
			a[L[k] - 1] = x;
			x = y;
			--L[k];
			++cnt;
		}
	}

	for (int j = 0; ++j < n;) {
		int value = a[j];
		i = j;
		while ((--i >= 0) && ((k = a[i]) > value))
			a[i + 1] = k;
		a[i + 1] = value;
	}
	delete[] L;
}

void Sort(int a[], int n, int sortType) {
	switch (sortType) {
	case 0: {
		SelectionSort(a, n);
		break;
	}
	case 1: {
		InsertionSort(a, n);
		break;
	}
	case 2: {
		BinaryInsertionSort(a, n);
		break;
	}
	case 3: {
		BubbleSort(a, n);
		break;
	}
	case 4: {
		ShakerSort(a, n);
		break;
	}
	case 5: {
		ShellSort(a, n);
		break;
	}
	case 6: {
		HeapSort(a, n);
		break;
	}
	case 7: {
		MergeSort(a, n);
		break;
	}
	case 8: {
		QuickSort(a, 0, n - 1);
		break;
	}
	case 9: {
		CountingSort(a, n);
		break;
	}
	case 10: {
		RadixSort(a, n);
		break;
	}
	case 11: {
		FlashSort(a, n);
		break;
	}
	}
}

