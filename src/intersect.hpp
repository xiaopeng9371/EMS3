#include <iostream>
#include <string>
#include<vector>
using namespace std;

#include <string.h>
#include <stdio.h>

//const string guard("TTTT");

void radix_sort(vector<string> & vCmbItem, vector<string> & vCmbTmp, int compact_count)
{
   // int compact_count = vCmbItem.size();
    int l = vCmbItem[0].size();
    for(int i = l - 1; i >= 0 ; i--)
    {
//        cout << "considering digit " << i << endl;
//        for (auto &s: vCmbItem)
//            cout << s << endl;
//        cout << "====================" << endl;

        int m1[4] = {0};
        for(unsigned int j = 0; j < compact_count;  j++ )
        {
            string s = vCmbItem[j];
            char x = s[i];

            switch (x) {
            case '*':
            {
                m1[0]++;
                m1[1]++;
                m1[2]++;
                m1[3]++;
                break;
            }
            case 'A':
            {
                m1[0]++;
                break;
            }
            case 'C':
            {
                m1[1]++;
                break;
            }
            case 'G':
            {
                m1[2]++;
                break;
            }
            case 'T':
            {
                m1[3]++;
                break;
            }
            default:
                break;
            }
        }
        m1[1] += m1[0];
        m1[2] += m1[1];
        m1[3] += m1[2];

        int expanded_count = m1[3];

//        vector<string> vCmbTmp;
//        vCmbTmp.resize(m1[3]);
        for(int j = compact_count-1; j >= 0; j--)
        {
            string s = vCmbItem[j];
            switch (s[i]) {
            case 'A':
                vCmbTmp[--m1[0]] = s;
                break;
            case 'C':
                vCmbTmp[--m1[1]] = s;
                break;
            case 'G':
                vCmbTmp[--m1[2]] = s;
                break;
            case 'T':
                vCmbTmp[--m1[3]] = s;
                break;
            case '*':
                s[i] = 'T'; vCmbTmp[--m1[3]] = s;
                s[i] = 'G'; vCmbTmp[--m1[2]] = s;
                s[i] = 'C'; vCmbTmp[--m1[1]] = s;
                s[i] = 'A'; vCmbTmp[--m1[0]] = s;
                break;
            default:
                break;
            }
        }
//        cout << "vcbtmp " << i << endl;
//        for (auto &s: vCmbTmp)
//            cout << s << endl;
//        cout << "--------------------" << endl;
        swap(vCmbItem, vCmbTmp);
        compact_count = expanded_count;
    }
}

//assume str1 is already sorted and without '*' and duplication free, str2 not necessarily sorted and may come with '*' and may have duplications
void intersect(vector<string> & str1,  vector<string> & str2, vector<string>& strTmp, int num_compact)
{
    str1.push_back(guard);

    radix_sort(str2, strTmp, num_compact);

//for (auto &s: str2)
//    cout << s << endl;

    auto first1 = str1.begin(), last1 = str1.end();
    auto first2 = str2.begin(), last2 = str2.end();
    auto result = first2, tmp = first2;
    while (first2 != last2)
    {
        while (*first1 < *first2) first1++;
        if (*first1 == *first2) {
            *result++ = std::move(*first2);
            tmp = first2++;
        }
        while (*tmp == *first2++);
    }
    str2.erase(result, last2);
}

/*
int main()
{
    vector<string> vStr1, vStr2, vStrTmp;
    const int ciLen1 = 4;
    const int ciLen2 = 6;
    const char* czValue[ciLen1] = {"AGT", "ATT", "TAG", "TCT"};
    for(int i = 0; i<ciLen1; i++)
    {
        vStr1.push_back(czValue[i]);
    }

    const char* czValue2[ciLen2] = {"AT*", "TCT", "*GT", "T*A", "TCT", "TTT"};
    for(int i = 0; i<ciLen2; i++)
    {
        vStr2.push_back(czValue2[i]);
    }

    int num_compact = vStr2.size();
    int num_expanded = 15;
    vStr2.resize(num_expanded);
    vStrTmp.resize(num_expanded);

    intersect(vStr1, vStr2, vStrTmp, num_compact);

    for(vector<string>::iterator itr = vStr2.begin(); itr != vStr2.end(); itr++)
    {
        cout << *itr << endl;
    }
}
*/
