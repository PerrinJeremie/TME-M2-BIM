%define name tv
%define version 0.5
%define release 1

Summary: TreeView X.
Name: %{name}
Version: %{version}
Release: %{release}
Source: http://darwin.zoology.gla.ac.uk/~rpage/treeviewx/%{name}-%{version}.tar.gz
Vendor: Roderic D. M. Page
URL: http://darwin.zoology.gla.ac.uk/~rpage/treeviewx/
License: GPL
Group: Sciences/Biology
BuildRoot: %{_builddir}/%{name}-buildroot
Prefix: %{_prefix}
Requires: wxGTK >= 2.4.0
BuildRequires: wxGTK-devel >= 2.4.0


%description
Phylogenetic tree viewer.


%prep
%setup -q

%build
if [ -x ./configure ]; then
  CFLAGS="$RPM_OPT_FLAGS" ./configure --prefix=%{_prefix}
else
  CFLAGS="$RPM_OPT_FLAGS" ./autogen.sh --prefix=%{_prefix}
fi
make

%install
rm -rf $RPM_BUILD_ROOT
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root)
%doc AUTHORS ChangeLog COPYING INSTALL NEWS README
%{_bindir}/*





%changelog
* Wed May 21 2003 R. D. M. Page
- First draft of the spec file
