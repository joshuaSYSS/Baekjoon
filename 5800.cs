using System;
using System.Linq;
using static System.Math;
public class baekjoon{
	public static void Main(){
		int n = Convert.ToInt32(Console.ReadLine());
		for(int i = 1;i <= n;i++){
		    string s = Console.ReadLine();
		    string[] sarr = s.Split(' ');
		    int[] arr = sarr
                        .Select(s2 => Convert.ToInt32(s2))
                        .ToArray();
            arr = arr.Skip(1).ToArray();
            Array.Sort(arr);
            int maxGap = 0;
            for(int j = 0;j < arr.Length - 1;j++){
                maxGap = Math.Max(maxGap, arr[j + 1] - arr[j]);
            }
            Console.WriteLine($"Class {i}");
            Console.WriteLine($"Max {arr[arr.Length - 1]}, Min {arr[0]}, Largest gap {maxGap}");
		}
	}
}
