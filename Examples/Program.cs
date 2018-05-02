using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

namespace Examples
{
    // Add this attribute to static methods to make
    // them available as examples.
    internal class ExampleMethodAttribute : Attribute { }

    class Program
    {
        private static MethodInfo[] GetExampleMethods () {
            Assembly assembly = Assembly.GetExecutingAssembly();
            MethodInfo[] methods = assembly.GetTypes()
                .SelectMany(t => t.GetMethods())
                .Where(m => m.GetCustomAttributes(typeof(ExampleMethodAttribute), false).Length > 0)
                .ToArray();
            return(methods);
        }


        static void Main(string[] args)
        {
            MethodInfo[] methods = GetExampleMethods();
            Dictionary<string, MethodInfo> index = new Dictionary<string, MethodInfo>();
            foreach (MethodInfo method in methods) {
                index.Add(method.Name, method);
            }

            if (args.Length == 0) {
                foreach(string key in index.Keys) {
                    Console.WriteLine(key);
                }
            } else {
                foreach (string arg in args) {
                    MethodInfo method = index[arg];
                    method.Invoke(null, null);
                }
            }

        }
    }
}
