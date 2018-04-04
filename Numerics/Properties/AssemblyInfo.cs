using System;
using System.Reflection;
using System.Resources;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Security;
#if !SILVERLIGHT
using System.Security.Permissions;
#endif

// General Information about an assembly is controlled through the following 
// set of attributes. Change these attribute values to modify the information
// associated with an assembly.
[assembly: AssemblyTitle("Meta.Numerics")]
[assembly: AssemblyDescription("")]
[assembly: AssemblyConfiguration("")]
[assembly: AssemblyCompany("Meta Numerics")]
[assembly: AssemblyCopyright("Copyright © Meta Numerics 2008-2018")]
[assembly: AssemblyTrademark("")]
[assembly: AssemblyCulture("")]

// CLS compliance
[assembly: CLSCompliant(true)]

// Define English as neutral resource language
// This attribute is required for the class to be used in Windows Universal Apps
[assembly: NeutralResourcesLanguage("en-us")]

// COM exposure
[assembly: ComVisible(false)]
[assembly: Guid("765978ee-7429-4132-a31f-233f8423f971")]

// Security
#if !SILVERLIGHT
//[assembly: SecurityPermission(SecurityAction.RequestMinimum, Execution = true)]
//[assembly: AllowPartiallyTrustedCallers]
#endif
//[assembly: SecurityTransparent]

// Version information for an assembly consists of the following four values:
//
//      Major Version
//      Minor Version 
//      Build Number
//      Revision
//
// You can specify all the values or you can default the Revision and Build Numbers 
// by using the '*' as shown below:
[assembly: AssemblyVersion("4.0.2.0")]
[assembly: AssemblyFileVersion("4.0.2.0")]

