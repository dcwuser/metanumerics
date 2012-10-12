using System;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Security;
using System.Security.Permissions;

// General Information about an assembly is controlled through the following 
// set of attributes. Change these attribute values to modify the information
// associated with an assembly.
[assembly: AssemblyTitle("Meta.Numerics")]
[assembly: AssemblyDescription("")]
[assembly: AssemblyConfiguration("")]
[assembly: AssemblyCompany("Meta Numerics")]
[assembly: AssemblyCopyright("Copyright © Meta Numerics 2008-2012")]
[assembly: AssemblyTrademark("")]
[assembly: AssemblyCulture("")]

// CLS compliance
[assembly: CLSCompliant(true)]

// COM exposure
[assembly: ComVisible(false)]
[assembly: Guid("765978ee-7429-4132-a31f-233f8423f971")]

// Security
#if !SILVERLIGHT
[assembly: SecurityPermission(SecurityAction.RequestMinimum, Execution = true)]
[assembly: AllowPartiallyTrustedCallers]
#endif

// Version information for an assembly consists of the following four values:
//
//      Major Version
//      Minor Version 
//      Build Number
//      Revision
//
// You can specify all the values or you can default the Revision and Build Numbers 
// by using the '*' as shown below:
[assembly: AssemblyVersion("2.2.0.0")]
[assembly: AssemblyFileVersion("2.2.0.0")]

